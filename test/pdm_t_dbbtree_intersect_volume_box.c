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
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_box_priv.h"
#include "pdm_box.h"
#include "pdm_dbbtree.h"

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const int    verbose = 1;
static const int    vtk     = 1;

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
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else
      _usage (EXIT_FAILURE);
    i++;
  }
}

/*
 * \brief Determine planes from box extents
 *
 * \param [in]  box_extents   Input box caracteristics
 * \param [out]  n             Table of normal vector (n = (a, b, c)) of each considered plane
 * \param [out]  plane_pt      Table of a point on plane (cartesian equation of plane being ax+by+cz+d=0 with d = -A.n) for each considered plane
 *
 */

static void
_box_extents_to_plane
(
double *box_extents,
double *n,
double *pt_plane
)
{
  double vect1[3];
  double vect2[3];
  double vect3[3];
  double vectn[3];
  double sign;

  double xmin = box_extents[0];
  double ymin = box_extents[1];
  double zmin = box_extents[2];
  double xmax = box_extents[3];
  double ymax = box_extents[4];
  double zmax = box_extents[5];

  // normal x

  vect1[0] = xmax - xmin;
  vect1[1] = 0;
  vect1[2] = 0;

  vect2[0] = 0;
  vect2[1] = 0;
  vect2[2] = zmax - zmin;

  PDM_CROSS_PRODUCT(vectn, vect1, vect2);

  vect3[0] = 0;
  vect3[1] = ymax - ymin;
  vect3[2] = 0;

  sign = PDM_SIGN(PDM_DOT_PRODUCT(vect3, vectn));

  n[0] = sign * vectn[0];
  n[1] = sign * vectn[1];
  n[2] = sign * vectn[2];

  pt_plane[0] = xmin;
  pt_plane[1] = ymin;
  pt_plane[2] = zmin;

  n[3] = - sign * vectn[0];
  n[4] = - sign * vectn[1];
  n[5] = - sign * vectn[2];

  pt_plane[3] = xmin;
  pt_plane[4] = ymax;
  pt_plane[5] = zmin;

  // normal y

  vect1[0] = 0;
  vect1[1] = ymax - ymin;
  vect1[2] = 0;

  vect2[0] = 0;
  vect2[1] = 0;
  vect2[2] = zmax - zmin;

  PDM_CROSS_PRODUCT(vectn, vect1, vect2);

  vect3[0] = xmax - xmin;
  vect3[1] = 0;
  vect3[2] = 0;

  sign = PDM_SIGN(PDM_DOT_PRODUCT(vect3, vectn));

  n[6] = sign * vectn[0];
  n[7] = sign * vectn[1];
  n[8] = sign * vectn[2];

  pt_plane[6] = xmin;
  pt_plane[7] = ymin;
  pt_plane[8] = zmin;

  n[9]  = - sign * vectn[0];
  n[10] = - sign * vectn[1];
  n[11] = - sign * vectn[2];

  pt_plane[9]  = xmax;
  pt_plane[10] = ymin;
  pt_plane[11] = zmin;

  // normal z

  vect1[0] = xmax - xmin;
  vect1[1] = 0;
  vect1[2] = 0;

  vect2[0] = 0;
  vect2[1] = ymax - ymin;
  vect2[2] = 0;

  PDM_CROSS_PRODUCT(vectn, vect1, vect2);

  vect3[0] = 0;
  vect3[1] = 0;
  vect3[2] = zmax - zmin;

  sign = PDM_SIGN(PDM_DOT_PRODUCT(vect3, vectn));

  n[12] = sign * vectn[0];
  n[13] = sign * vectn[1];
  n[14] = sign * vectn[2];

  pt_plane[12] = xmin;
  pt_plane[13] = ymin;
  pt_plane[14] = zmin;

  n[15] = - sign * vectn[0];
  n[16] = - sign * vectn[1];
  n[17] = - sign * vectn[2];

  pt_plane[15] = xmin;
  pt_plane[16] = ymin;
  pt_plane[17] = zmax;

}

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  int n_part = 1;

  _read_args (argc,
              argv);

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  // Create dbbtree boxes

  int          n_dbbtree_boxes     = 2;
  double      *dbbtree_box_extents = malloc(sizeof(double)      * n_dbbtree_boxes * 6);
  PDM_g_num_t *dbbtree_box_g_num   = malloc(sizeof(PDM_g_num_t) * n_dbbtree_boxes);

  for (int j = 0; j < n_dbbtree_boxes; j++) {
    dbbtree_box_extents[0 + 6 * j] = 0   + j * 0.5 + i_rank * 1;
    dbbtree_box_extents[1 + 6 * j] = 0;
    dbbtree_box_extents[2 + 6 * j] = 0;
    dbbtree_box_extents[3 + 6 * j] = 0.5 + j * 0.5 + i_rank * 1;
    dbbtree_box_extents[4 + 6 * j] = 1;
    dbbtree_box_extents[5 + 6 * j] = 1;

    dbbtree_box_g_num[j] = (i_rank * n_dbbtree_boxes) + j + 1;
  }

  if (vtk) {
    char filename1[999];
    sprintf(filename1, "dbbtree_boxes_%d.vtk", i_rank);
    PDM_vtk_write_boxes(filename1,
                        n_dbbtree_boxes,
                        dbbtree_box_extents,
                        dbbtree_box_g_num);
  }

  // Create volume boxes

  int          n_volume_boxes     = 1;
  double      *volume_box_extents = malloc(sizeof(double)      * n_volume_boxes * 6);
  PDM_g_num_t *volume_box_g_num   = malloc(sizeof(PDM_g_num_t) * n_volume_boxes);

  for (int j = 0; j < n_volume_boxes; j++) {
    volume_box_extents[0 + 6 * j] = 0.5 + j * 0.5 + i_rank * 1;
    volume_box_extents[1 + 6 * j] = 0;
    volume_box_extents[2 + 6 * j] = 0;
    volume_box_extents[3 + 6 * j] = 0.5 + j * 0.5 + (i_rank+1) * 1;
    volume_box_extents[4 + 6 * j] = 1;
    volume_box_extents[5 + 6 * j] = 1;

    volume_box_g_num[j] = (i_rank * n_volume_boxes) + j + 1;
  }

  if (vtk) {
    char filename2[999];
    sprintf(filename2, "volume_boxes_%d.vtk", i_rank);
    PDM_vtk_write_boxes(filename2,
                        n_volume_boxes,
                        volume_box_extents,
                        volume_box_g_num);
  }

  // Create dbbtree

  const int dim = 3;
  double l_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL,
                         -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int i = 0; i < n_dbbtree_boxes; i++) {
    for (int k = 0; k < 3; k++) {
      l_extents[k]     = PDM_MIN(l_extents[k],     dbbtree_box_extents[6*i + k]);
      l_extents[k + 3] = PDM_MAX(l_extents[k + 3], dbbtree_box_extents[6*i + k + 3]);
    }
  }

  double g_extents[6];
  PDM_MPI_Allreduce(l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX(max_range, g_extents[i+3] - g_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    g_extents[i]   -= max_range * 1.1e-3;
    g_extents[i+3] += max_range * 1.0e-3;
  }

  PDM_dbbtree_t *dbbt = PDM_dbbtree_create(comm, dim, g_extents);

  PDM_box_set_t *box_set = PDM_dbbtree_boxes_set(dbbt,
                                                 n_part,
                                                 &n_dbbtree_boxes,
                                                 (const double **) &dbbtree_box_extents,
                                                 (const PDM_g_num_t **) &dbbtree_box_g_num);

  int *volume_plane_idx = PDM_array_new_idx_from_const_stride_int(6, n_volume_boxes);

  double *plane_normal   = malloc(sizeof(double) * 3 * n_volume_boxes * 6);
  double *plane_pt_coord = malloc(sizeof(double) * 3 * n_volume_boxes * 6);

  for (int ibox = 0; ibox < n_volume_boxes; ibox++) {
  _box_extents_to_plane(volume_box_extents + 6 * ibox,
                        plane_normal       + 6 * ibox * 3,
                        plane_pt_coord     + 6 * ibox * 3);
  } // end loop on volume boxes

  int         *volume_boxes_idx   = NULL;
  PDM_g_num_t *volume_boxes_g_num = NULL;
  PDM_dbbtree_volumes_intersect_boxes(dbbt,
                                      n_volume_boxes,
                                      volume_box_g_num,
                                      volume_plane_idx,
                                      plane_normal,
                                      plane_pt_coord,
                                      &volume_boxes_idx,
                                      &volume_boxes_g_num);

  if (verbose) {
    log_trace("VOLUME-BOX INTERSECTION\n");
    for (int ivol = 0; ivol < n_volume_boxes; ivol++) {
      log_trace("--> volume "PDM_FMT_G_NUM" is intersected by ");
      for (int i = volume_boxes_idx[ivol]; i < volume_boxes_idx[ivol+1]; i++) {
        log_trace("%d ", volume_boxes_g_num[i]);
      }
      log_trace("\n");
    }
  }

  PDM_dbbtree_free(dbbt);
  PDM_box_set_destroy(&box_set);

  PDM_MPI_Finalize ();
}

