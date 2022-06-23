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
  double direction_pt[3] = {1, 20, 0}; // C---D-->
  double theta = PDM_PI / 3;
  double eps   = 1;

  double box_extents[6] = {-1, 1, -1, 3, 5, 3}; // xmin, ymin, zmin, xmax, ymax, zmax
  double n[12];
  double pt_plane[12];

  // Determine eps translation planes
  // B--{eps}--G
  double CB[3] = {edge[3]-edge[6], edge[4]-edge[7], edge[5]-edge[8]};
  double inverse_module_CB = 1 / PDM_MODULE(CB);
  pt_plane[0] = edge[3] + (1+eps) * CB[0]; // * inverse_module_CB;
  pt_plane[1] = edge[4] + (1+eps) * CB[1]; // * inverse_module_CB;
  pt_plane[2] = edge[5] + (1+eps) * CB[2]; // * inverse_module_CB;
  n[0] = -CB[0];
  n[1] = -CB[1];
  n[2] = -CB[2];

  // A--{eps}--H
  double CA[3] = {edge[3]-edge[0], edge[4]-edge[1], edge[5]-edge[2]};
  double inverse_module_CA = 1 / PDM_MODULE(CA);
  pt_plane[3] = edge[3] + (1+eps) * CA[0]; // * inverse_module_CA;
  pt_plane[4] = edge[4] + (1+eps) * CA[1]; // * inverse_module_CA;
  pt_plane[5] = edge[5] + (1+eps) * CA[2]; // * inverse_module_CA;
  n[3] = -CA[0];
  n[4] = -CA[1];
  n[5] = -CA[2];

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
  double CE[3] = {pt_plane[6]-edge[3], pt_plane[7]-edge[4], pt_plane[8]-edge[5]};
  PDM_CROSS_PRODUCT(prod_vect, CA, CE);
  double ED[3] = {direction_pt[0] - pt_plane[6], direction_pt[1] - pt_plane[7], direction_pt[2] - pt_plane[8]};
  double sign = PDM_SIGN(PDM_DOT_PRODUCT(prod_vect, ED));
  n[6] = sign * prod_vect[0];
  n[7] = sign * prod_vect[1];
  n[8] = sign * prod_vect[2];

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
  double CF[3] = {pt_plane[9]-edge[4], pt_plane[10]-edge[4], pt_plane[11]-edge[5]};
  PDM_CROSS_PRODUCT(prod_vect, CA, CF);
  double FD[3] = {direction_pt[0] - pt_plane[9], direction_pt[1] - pt_plane[10], direction_pt[2] -  pt_plane[11]};
  sign = PDM_SIGN(PDM_DOT_PRODUCT(prod_vect, FD));
  n[9]  = sign * prod_vect[0];
  n[10] = sign * prod_vect[1];
  n[11] = sign * prod_vect[2];

  // Check if box is in volume
  int check = _box_in_volume(4, n, pt_plane, box_extents);
  log_trace("box is in volume = %d\n", check);

  // vtk output of atomic test case

  char *filename1 = "box.vtk";
  PDM_g_num_t *box_g_num = malloc(sizeof(PDM_g_num_t) * 1);
  box_g_num[0] = 1;

  PDM_vtk_write_boxes(filename1,
                      1,
                      box_extents,
                      box_g_num);

  char *filename2 = "line.vtk";
  double *coord = malloc(sizeof(double) * 6);
  coord[0] = edge[3];
  coord[1] = edge[4];
  coord[2] = edge[5];
  coord[3] = direction_pt[0];
  coord[4] = direction_pt[1];
  coord[5] = direction_pt[2];
  PDM_g_num_t *line_g_num = malloc(sizeof(PDM_g_num_t) * 1);
  line_g_num[0] = 1;

  PDM_vtk_write_lines(filename2,
                      1,
                      coord,
                      line_g_num,
                      NULL);

  char *filename3 = "planes.vtk";
  double *vtx_coord = malloc(sizeof(double) * 30);
  PDM_g_num_t *vtx_g_num = malloc(sizeof(PDM_g_num_t) * 10);
  int *face_vtx = malloc(sizeof(int) * 12);

  // A
  vtx_coord[0] = edge[0];
  vtx_coord[1] = edge[1];
  vtx_coord[2] = edge[2];

  // B
  vtx_coord[3] = edge[6];
  vtx_coord[4] = edge[7];
  vtx_coord[5] = edge[8];

  // E
  vtx_coord[6] = pt_plane[6];
  vtx_coord[7] = pt_plane[7];
  vtx_coord[8] = pt_plane[8];

  // F
  vtx_coord[9]  = pt_plane[9];
  vtx_coord[10] = pt_plane[10];
  vtx_coord[11] = pt_plane[11];

  // G
  vtx_coord[12] = pt_plane[0];
  vtx_coord[13] = pt_plane[1];
  vtx_coord[14] = pt_plane[2];

  // H
  vtx_coord[15] = pt_plane[3];
  vtx_coord[16] = pt_plane[4];
  vtx_coord[17] = pt_plane[5];

  double CH[3] = {pt_plane[3]-edge[3], pt_plane[4]-edge[4], pt_plane[5]-edge[5]};
  double module_CH = PDM_MODULE(CH);

  // I
  vtx_coord[18] = direction_pt[0] + module_CH * CA[0] * inverse_module_CA;
  vtx_coord[19] = direction_pt[1] + module_CH * CA[1] * inverse_module_CA;
  vtx_coord[20] = direction_pt[2] + module_CH * CA[2] * inverse_module_CA;

  // J
  vtx_coord[21] = pt_plane[6] + module_CH * CA[0] * inverse_module_CA;
  vtx_coord[22] = pt_plane[7] + module_CH * CA[1] * inverse_module_CA;
  vtx_coord[23] = pt_plane[8] + module_CH * CA[2] * inverse_module_CA;

  double CG[3] = {pt_plane[0]-edge[3], pt_plane[1]-edge[4], pt_plane[2]-edge[5]};
  double module_CG = PDM_MODULE(CG);

  // K
  vtx_coord[24] = direction_pt[0] + module_CG * CB[0] * inverse_module_CB;
  vtx_coord[25] = direction_pt[1] + module_CG * CB[1] * inverse_module_CB;
  vtx_coord[26] = direction_pt[2] + module_CG * CB[2] * inverse_module_CB;

  // L
  vtx_coord[27] = pt_plane[6] + module_CG * CB[0] * inverse_module_CB;
  vtx_coord[28] = pt_plane[7] + module_CG * CB[1] * inverse_module_CB;
  vtx_coord[29] = pt_plane[8] + module_CG * CB[2] * inverse_module_CB;

  for (int i = 0; i < 10; i++) {
    vtx_g_num[i] = i + 1;
  }

  face_vtx[0]  = 1;
  face_vtx[1]  = 2;
  face_vtx[2]  = 3;
  face_vtx[3]  = 1;
  face_vtx[4]  = 2;
  face_vtx[5]  = 4;
  face_vtx[6]  = 5;
  face_vtx[7]  = 9;
  face_vtx[8]  = 10;
  face_vtx[9]  = 6;
  face_vtx[10] = 7;
  face_vtx[11] = 8;


  PDM_vtk_write_std_elements(filename3,
                             10,
                             vtx_coord,
                             vtx_g_num,
                             PDM_MESH_NODAL_TRIA3,
                             4,
                             face_vtx,
                             NULL,
                             0,
                             NULL,
                             NULL);

  char *filename4 = "normal.vtk";

  double *vector_normal[1] = {n};

  const char* normal_name[] = {"n", 0};

  PDM_vtk_write_point_cloud_with_field(filename4,
                                       4,
                                       pt_plane,
                                       NULL,
                                       NULL,
                                       0,
                                       NULL,
                                       NULL,
                                       1,
                       (const char **) &normal_name,
                     (const double **) &vector_normal,
                                       0,
                                       NULL,
                                       NULL);

  /* dbbtree test case */

  PDM_MPI_Finalize ();

}
