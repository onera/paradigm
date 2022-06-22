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
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"

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
 * \brief Determine if a box is on the side of the plane where the normal points to
 *
 * \param [in]  n                Normal vector (n = (a, b, c))
 * \param [in]  plane_pt         Point on plane (cartesian equation of plane being ax+by+cz+d=0 with d = -A.n)
 * \param [in]  box_extents      Box to determine
 *
 * \return 1 if the box is on the side of the plane where the normal points to, 0 otherwise
 */

static int
_plane_box_side
(
double *n,
double *plane_pt,
double *box_extents
)
{
  double box_pt[3];
  double vect[3];

  for (int x = 0; x < 4; x += 3) {
    box_pt[0] = box_extents[x];
    for (int y = 1; y < 5; y += 3) {
      box_pt[1] = box_extents[y];
      for (int z = 2; z < 6; z += 3) {
        box_pt[2] = box_extents[z];
        vect[0] = box_pt[0] - plane_pt[0]; vect[1] = box_pt[1] - plane_pt[1]; vect[2] = box_pt[2] - plane_pt[2];
        if (PDM_DOT_PRODUCT(vect, n) > 0) { // if >= 0 also considers when point is on the plane
          return 1;
        }
      } // end loop on z
    } // end loop on y
  } // end loop on x
  return 0;
}

/*
 * \brief Determine a box is in a given volume region
 *
 * \param [in]  n_planes      Number of planes difining the volume
 * \param [in]  n             Table of normal vector (n = (a, b, c)) of each considered plane
 * \param [in]  plane_pt      Table of a point on plane (cartesian equation of plane being ax+by+cz+d=0 with d = -A.n) for each considered plane
 * \param [in]  box_extents   Point to determine
 *
 * \return 1 if the plane_pt is on the side of the plane where the normal points to, 0 otherwise
 */


static int
_box_in_volume
(
int      n_planes,
double  *n,
double  *plane_pt,
double  *box_extents
)
{
  double n_iplane[3];
  double plane_pt_iplane[3];

  for (int iplane = 0; iplane < n_planes; iplane++) {

    n_iplane[0] = n[3*iplane];
    n_iplane[1] = n[3*iplane+1];
    n_iplane[2] = n[3*iplane+2];

    plane_pt_iplane[0] = plane_pt[3*iplane];
    plane_pt_iplane[1] = plane_pt[3*iplane+1];
    plane_pt_iplane[2] = plane_pt[3*iplane+2];

    if (_plane_box_side(n_iplane, plane_pt_iplane, box_extents) == 0) {
      return 0;
    }
  } // end loop on planes
  return 1;
}

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Init (&argc, &argv);

  int           i_rank;
  int           numProcs;

  _read_args (argc,
              argv);

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  /* TO DO :
   *   - create a PDM_box_tree_inside_volume_boxes
   *   - create _volume_intersect_shared_box_tree
   *   - create PDM_dbbtree_boxes_inside_volume
   */

  /* Atomic test case */

  // Set up

  double edge[9] = {0, 0, 0, 1, 0, 0, 2, 0, 0}; // A---C---B
  double direction_pt[3] = {1, 1, 0}; // C---D-->
  double theta = PDM_PI / 3;
  double eps   = 1;

  double box_extents[6] = {-1, 1, -1, 3, 5, 3}; // xmin, ymin, zmin, xmax, ymax, zmax
  double n[12];
  double pt_plane[12];

  // Determine eps translation planes
  // B--{eps}--G
  double BC[3] = {edge[6]-edge[3], edge[7]-edge[4], edge[8]-edge[5]};
  double inverse_module_BC = 1 / PDM_MODULE(BC);
  pt_plane[0] = edge[6] + eps * BC[0] * inverse_module_BC;
  pt_plane[1] = edge[7] + eps * BC[1] * inverse_module_BC;
  pt_plane[2] = edge[8] + eps * BC[2] * inverse_module_BC;
  n[0] = BC[0];
  n[1] = BC[1];
  n[2] = BC[2];

  // A--{eps}--H
  double CA[3] = {edge[3]-edge[0], edge[4]-edge[1], edge[5]-edge[2]};
  double inverse_module_CA = 1 / PDM_MODULE(CA);
  pt_plane[3] = edge[6] + eps * CA[0] * inverse_module_CA;
  pt_plane[4] = edge[7] + eps * CA[1] * inverse_module_CA;
  pt_plane[5] = edge[8] + eps * CA[2] * inverse_module_CA;
  n[3] = CA[0];
  n[4] = CA[1];
  n[5] = CA[2];

  // Determine theta angle planes E---D---F
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);


  double AB[3] = {edge[6]-edge[0], edge[7]-edge[1], edge[8]-edge[2]};
  double inverse_module_AB = 1 / PDM_MODULE(AB);
  double AB_normalised[3] = {AB[0] * inverse_module_AB, AB[1] * inverse_module_AB, AB[2] * inverse_module_AB};

  pt_plane[6]  = (cos_theta + (1 - cos_theta) * AB_normalised[0] * AB_normalised[0]) * direction_pt[0];
  pt_plane[6] += (AB_normalised[0] * AB_normalised[1] * (1 - cos_theta) - AB_normalised[2] * sin_theta) * direction_pt[1];
  pt_plane[6] += (AB_normalised[0] * AB_normalised[2] * (1 - cos_theta) + AB_normalised[1] * sin_theta) * direction_pt[2];
  pt_plane[7]  = (AB_normalised[1] * AB_normalised[0] * (1 - cos_theta) + AB_normalised[2] * sin_theta) * direction_pt[0];
  pt_plane[7] += (cos_theta + (1 - cos_theta) * AB_normalised[1] * AB_normalised[1]) * direction_pt[1];
  pt_plane[7] += (AB_normalised[1] * AB_normalised[2] * (1 - cos_theta) - AB_normalised[0] * sin_theta) * direction_pt[2];
  pt_plane[8]  = (AB_normalised[2] * AB_normalised[0] * (1 - cos_theta) - AB_normalised[1] * sin_theta) * direction_pt[0];
  pt_plane[8] += (AB_normalised[2] * AB_normalised[1] * (1 - cos_theta) + AB_normalised[0] * sin_theta) * direction_pt[1];
  pt_plane[8] += (cos_theta + (1 - cos_theta) * AB_normalised[2] * AB_normalised[2]) * direction_pt[2];
  double prod_vect[3];
  double CE[3] = {pt_plane[6]-edge[0], pt_plane[7]-edge[1], pt_plane[8]-edge[2]};
  PDM_CROSS_PRODUCT(prod_vect, CA, CE);
  double ED[3] = {pt_plane[6] - direction_pt[0], pt_plane[7] - direction_pt[1], pt_plane[8] - direction_pt[2]};
  n[6] = PDM_SIGN(ED[0]) * prod_vect[0];
  n[7] = PDM_SIGN(ED[1]) * prod_vect[1];
  n[8] = PDM_SIGN(ED[2]) * prod_vect[2];

  double cos_minus_theta = cos(-theta);
  double sin_minus_theta = sin(-theta);

  pt_plane[9]   = (cos_minus_theta + (1 - cos_minus_theta) * AB_normalised[0] * AB_normalised[0]) * direction_pt[0];
  pt_plane[9]  += (AB_normalised[0] * AB_normalised[1] * (1 - cos_minus_theta) - AB_normalised[2] * sin_minus_theta) * direction_pt[1];
  pt_plane[9]  += (AB_normalised[0] * AB_normalised[2] * (1 - cos_minus_theta) + AB_normalised[1] * sin_minus_theta) * direction_pt[2];
  pt_plane[10]  = (AB_normalised[1] * AB_normalised[0] * (1 - cos_minus_theta) + AB_normalised[2] * sin_minus_theta) * direction_pt[0];
  pt_plane[10] += (cos_minus_theta + (1 - cos_minus_theta) * AB_normalised[1] * AB_normalised[1]) * direction_pt[1];
  pt_plane[10] += (AB_normalised[1] * AB_normalised[2] * (1 - cos_minus_theta) - AB_normalised[0] * sin_minus_theta) * direction_pt[2];
  pt_plane[11]  = (AB_normalised[2] * AB_normalised[0] * (1 - cos_minus_theta) - AB_normalised[1] * sin_minus_theta) * direction_pt[0];
  pt_plane[11] += (AB_normalised[2] * AB_normalised[1] * (1 - cos_minus_theta) + AB_normalised[0] * sin_minus_theta) * direction_pt[1];
  pt_plane[11] += (cos_minus_theta + (1 - cos_minus_theta) * AB_normalised[2] * AB_normalised[2]) * direction_pt[2];
  n[9]  = -n[6];
  n[10] = -n[7];
  n[11] = -n[8];

  // Check if box is in volume
  int check = _box_in_volume(4, n, pt_plane, box_extents);
  log_trace("box is in volume = %d\n", check);

  /* dbbtree test case */

  PDM_MPI_Finalize ();

}
