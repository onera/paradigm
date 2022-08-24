#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_dbbtree.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_multipart.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_plane.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_triangle.h"
#include "pdm_distrib.h"
#include "pdm_unique.h"

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const int verbose = 0;
static const int vtk     = 1;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -h               This message.\n\n");


  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *n_back,
 PDM_g_num_t   *n_work
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-nb") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long n = atol(argv[i]);
        *n_back = (PDM_g_num_t) n;
      }
    }

    else if (strcmp(argv[i], "-nw") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long n = atol(argv[i]);
        *n_work = (PDM_g_num_t) n;
      }
    }

    else
      _usage (EXIT_FAILURE);
    i++;
  }
}





/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  int n_part = 1;
  PDM_g_num_t n_back = 10;
  PDM_g_num_t n_work = 1;
 _read_args (argc,
              argv,
              &n_back,
              &n_work);

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  const double x_center = 0.1;
  const double y_center = 0.2;
  const double z_center = 0.3;
  const double radius   = 3.14;


  /*
   *  Generate a fine background distributed mesh (single surface group)
   */
  double      *dback_vtx_coord    = NULL;
  int         *dback_face_vtx_idx = NULL;
  PDM_g_num_t *dback_face_vtx     = NULL;
  PDM_g_num_t *back_distrib_vtx   = NULL;
  PDM_g_num_t *back_distrib_face  = NULL;

  // PDM_sphere_surf_gen(comm,
  //                     2*n_back,
  //                     n_back,
  //                     x_center,
  //                     y_center,
  //                     z_center,
  //                     radius,
  //                     &dback_vtx_coord,
  //                     &dback_face_vtx_idx,
  //                     &dback_face_vtx,
  //                     &back_distrib_vtx,
  //                     &back_distrib_face);

  PDM_sphere_surf_icosphere_gen(comm,
                                n_back,
                                x_center,
                                y_center,
                                z_center,
                                radius,
                                &dback_vtx_coord,
                                &dback_face_vtx_idx,
                                &dback_face_vtx,
                                &back_distrib_vtx,
                                &back_distrib_face);

  /*
   *  Build the back mesh dbbtree
   */

  /* Assemble partitions from block-distribution of faces */
  int dback_n_face = back_distrib_face[i_rank+1] - back_distrib_face[i_rank];

  PDM_g_num_t *dback_face_ln_to_gn = malloc(dback_n_face * sizeof(PDM_g_num_t));
  for (int i = 0; i < dback_n_face; ++i) {
    dback_face_ln_to_gn[i] = back_distrib_face[i_rank] + i + 1;
  }


  PDM_g_num_t *pback_vtx_ln_to_gn = NULL;
  int         *pback_face_vtx_idx = NULL;
  int         *pback_face_vtx     = NULL;
  int          pback_n_vtx;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           back_distrib_face,
                                                           dback_face_vtx_idx,
                                                           dback_face_vtx,
                                                           dback_n_face,
                                  (const PDM_g_num_t *)    dback_face_ln_to_gn,
                                                           &pback_n_vtx,
                                                           &pback_vtx_ln_to_gn,
                                                           &pback_face_vtx_idx,
                                                           &pback_face_vtx);

  double **tmp_pback_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        n_part,
                                        back_distrib_vtx,
                                        dback_vtx_coord,
                                        &pback_n_vtx,
                 (const PDM_g_num_t **) &pback_vtx_ln_to_gn,
                                        &tmp_pback_vtx_coord);
  double *pback_vtx_coord = tmp_pback_vtx_coord[0];
  free(tmp_pback_vtx_coord);



  /* Compute the bounding boxes of local faces */
  double *back_face_extents = malloc(sizeof(double) * dback_n_face * 6);
  const double eps_extents = 1.0e-6;
  for (int iface = 0; iface < dback_n_face; iface++) {

    double *tmp_extents = back_face_extents + 6*iface;
    for (int k = 0; k < 3; k++) {
      tmp_extents[k]     =  HUGE_VAL;
      tmp_extents[3 + k] = -HUGE_VAL;
    }

    for (int ivtx = pback_face_vtx_idx[iface]; ivtx < pback_face_vtx_idx[iface+1]; ivtx++) {
      int vtx_id = pback_face_vtx[ivtx]-1;
      double *tmp_coord = pback_vtx_coord + 3*vtx_id;
      for (int k = 0; k < 3; k++) {
        tmp_extents[k]     =  PDM_MIN(tmp_coord[k], tmp_extents[k]);
        tmp_extents[3 + k] =  PDM_MAX(tmp_coord[k], tmp_extents[3 + k]);
      }
    } // end loop on vertices of iface

    // Threshold/Inflate extents
    for (int k = 0; k < 3; k++) {
      double size = tmp_extents[3 + k] - tmp_extents[k];

      // tmp_extents[k]     -= size*eps_extents;
      // tmp_extents[3 + k] += size*eps_extents;

      if (size < eps_extents) {
        tmp_extents[k]     -= 0.5*eps_extents;
        tmp_extents[3 + k] += 0.5*eps_extents;
      }
    }

  } // end loop on background faces


  if (vtk) {
    char filename[999];

    sprintf(filename, "back_faces_%d.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           pback_n_vtx,
                           pback_vtx_coord,
                           pback_vtx_ln_to_gn,
                           dback_n_face,
                           pback_face_vtx_idx,
                           pback_face_vtx,
                           dback_face_ln_to_gn,
                           NULL);

    sprintf(filename, "back_face_extents_%d.vtk", i_rank);
    PDM_vtk_write_boxes(filename,
                        dback_n_face,
                        back_face_extents,
                        dback_face_ln_to_gn);
  }
  free(pback_vtx_ln_to_gn);
  free(pback_vtx_coord);
  free(pback_face_vtx_idx);
  free(pback_face_vtx);


  /* Create and build dbbtree */
  const int dim = 3;
  double l_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL,
                         -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int i = 0; i < dback_n_face; i++) {
    for (int k = 0; k < 3; k++) {
      l_extents[k]     = PDM_MIN(l_extents[k],     back_face_extents[6*i + k]);
      l_extents[k + 3] = PDM_MAX(l_extents[k + 3], back_face_extents[6*i + k + 3]);
    }
  }

  double g_extents[6];
  PDM_MPI_Allreduce(l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX(max_range, g_extents[i+3] - g_extents[i]);
  }
  // Inflate and break symmetry
  for (int i = 0; i < 3; i++) {
    g_extents[i]   -= max_range * 1.1e-3;
    g_extents[i+3] += max_range * 1.0e-3;
  }

  PDM_dbbtree_t *dbbt = PDM_dbbtree_create(comm, dim, g_extents);

  PDM_box_set_t *box_set = PDM_dbbtree_boxes_set(dbbt,
                                                 n_part,
                                                 &dback_n_face,
                               (const double **) &back_face_extents,
                          (const PDM_g_num_t **) &dback_face_ln_to_gn);
  free(back_face_extents);
  free(dback_face_ln_to_gn);




  /*
   *  Generate a coarse 'work' partitioned mesh
   */

  /* Generate block-distributed nodal mesh */
  PDM_dmesh_nodal_t *dmn = NULL;

  // PDM_sphere_surf_gen_nodal(comm,
  //                           2*n_work,
  //                           n_work,
  //                           x_center,
  //                           y_center,
  //                           z_center,
  //                           radius,
  //                           &dmn);
  PDM_sphere_surf_icosphere_gen_nodal(comm,
                                      n_work,
                                      x_center,
                                      y_center,
                                      z_center,
                                      radius,
                                      &dmn);

  /* Split the mesh */
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
  int n_zone                   = 1;
  int *n_part_zones            = (int *) malloc(sizeof(int) * n_zone);
  n_part_zones[0]              = n_part;

  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                n_part_zones,
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

  free(n_part_zones);



  // Vertices
  int     i_part          = 0;
  double *pwork_vtx_coord = NULL;
  int pwork_n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                     0,
                                                     i_part,
                                                     &pwork_vtx_coord,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_g_num_t *pwork_vtx_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_VERTEX,
                                  &pwork_vtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  // Edges
  int *pwork_edge_vtx     = NULL;
  int *pwork_edge_vtx_idx = NULL;
  int  pwork_n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                          0,
                                                          i_part,
                                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                          &pwork_edge_vtx,
                                                          &pwork_edge_vtx_idx,
                                                          PDM_OWNERSHIP_KEEP);

  if (pwork_edge_vtx_idx != NULL) free(pwork_edge_vtx_idx);

  PDM_g_num_t *pwork_edge_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_EDGE,
                                  &pwork_edge_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);


  if (vtk) {
    char filename[999];

    sprintf(filename, "work_edges_%d.vtk", i_rank);
    PDM_vtk_write_std_elements(filename,
                               pwork_n_vtx,
                               pwork_vtx_coord,
                               pwork_vtx_ln_to_gn,
                               PDM_MESH_NODAL_BAR2,
                               pwork_n_edge,
                               pwork_edge_vtx,
                               pwork_edge_ln_to_gn,
                               0,
                               NULL,
                               NULL);
  }



  /*
   *  Free memory
   */
  PDM_multipart_free(mpart);
  PDM_DMesh_nodal_free(dmn);
  free(dback_vtx_coord);
  free(dback_face_vtx_idx);
  free(dback_face_vtx);
  free(back_distrib_vtx);
  free(back_distrib_face);

  PDM_dbbtree_free(dbbt);
  PDM_box_set_destroy(&box_set);


  PDM_MPI_Finalize ();

  return 0;
}
