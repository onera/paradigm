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
#include "pdm_ho_bezier_basis.h"
#include "pdm_triangle.h"

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

  // Set-up test
  // double *line       = malloc(sizeof(double) * 6);
  double tria_coord[21];
  // double ip[3];

  // line[0] = 0; line[3] = 0;
  // line[1] = 0; line[4] = 0;
  // line[2] = 0; line[5] = 0;

  tria_coord[0] = 0; tria_coord[3] = 20; tria_coord[6] = 0;
  tria_coord[1] = 0; tria_coord[4] = 0; tria_coord[7] = 20;
  tria_coord[2] = 0; tria_coord[5] = 0; tria_coord[8] = 0;

  tria_coord[9]  = 10; tria_coord[12] = 0;  tria_coord[15] = 10;
  tria_coord[10] = 0;  tria_coord[13] = 10; tria_coord[16] = 10;
  tria_coord[11] = 0;  tria_coord[14] = 0;  tria_coord[17] = 0;

  // // Do triangle P1 line intersection

  // PDM_triangle_line_intersection(line, tria_coord, ip);

  // // Construct matrices.  Since we have over determined system, need to find
  // // which 2 out of 3 equations to use to develop equations. (Any 2 should
  // // work since we've projected point to plane.)

  // double rhs[2], c1[2], c2[2];

  // for (int i = 0; i < 2; i++) {
  //   rhs[i] = ip[i] - tria_coord[6 + i];
  //   c1[i] = tria_coord[i] - tria_coord[6 + i];
  //   c2[i] = tria_coord[3 + i] - tria_coord[6 + i];
  // }

  // double det = PDM_DETERMINANT2X2(c1,c2);

  // double pcoords[3];

  // pcoords[0] = PDM_DETERMINANT2X2(rhs,c2) / det;
  // pcoords[1] = PDM_DETERMINANT2X2(c1,rhs) / det;

  // double weights[3];
  // double closest_point[3];

  // weights[0] = 1 - (pcoords[0] + pcoords[1]);
  // weights[1] = pcoords[0];
  // weights[2] = pcoords[1];
  // if ( weights[0] >= 0.0 && weights[0] <= 1.0 &&
  //      weights[1] >= 0.0 && weights[1] <= 1.0 &&
  //      weights[2] >= 0.0 && weights[2] <= 1.0 ) {

  //   // Projection distance

  //   closest_point[0] = ip[0];
  //   closest_point[1] = ip[1];
  //   closest_point[2] = ip[2];

  // } // end is inside case

  // // Get uvw in triangle Pn

  // double uv_clossest_Pn[3];

  // double _uvPn_sub_tria[6];

  // _uvPn_sub_tria[0] = uv_nodes[2*idx1];
  // _uvPn_sub_tria[1] = uv_nodes[2*idx1+1];
  // _uvPn_sub_tria[2] = uv_nodes[2*idx2];
  // _uvPn_sub_tria[3] = uv_nodes[2*idx2+1];
  // _uvPn_sub_tria[4] = uv_nodes[2*idx3];
  // _uvPn_sub_tria[5] = uv_nodes[2*idx3+1];

  // for (int j = 0; j < 2; j++) {
  //   for (int k = 0; k < 3; k++) {
  //     uv_clossest_Pn[j] += weights[k] * uvPn_sub_tria[2*k + j];
  //   }
  // }

  // Get xyz coordinates from uvw coordinates

  double uvw[3];
  double weights[6];

  uvw[0] = 0.25;
  uvw[1] = 0.25;
  uvw[2] = 0;

  PDM_ho_bezier_basis(PDM_MESH_NODAL_TRIAHO,
                      2,
                      1,
                      uvw,
                      weights);

  double xyz_Pn[3] = {0, 0, 0};

  for (int j = 0; j < 3; j++) {
    xyz_Pn[j] += weights[0] * tria_coord[3*0+j]; // (i,j,k)=(0,0,2)
    xyz_Pn[j] += weights[1] * tria_coord[3*3+j]; // (i,j,k)=(1,0,1)
    xyz_Pn[j] += weights[2] * tria_coord[3*1+j]; // (i,j,k)=(2,0,0)
    xyz_Pn[j] += weights[3] * tria_coord[3*4+j]; // (i,j,k)=(0,1,1)
    xyz_Pn[j] += weights[4] * tria_coord[3*5+j]; // (i,j,k)=(1,1,0)
    xyz_Pn[j] += weights[5] * tria_coord[3*2+j]; // (i,j,k)=(0,2,0)
  }

  PDM_log_trace_array_double(xyz_Pn, 3, "xyz on P2: ");

  // vtk output

  tria_coord[18] = xyz_Pn[0];
  tria_coord[19] = xyz_Pn[1];
  tria_coord[20] = xyz_Pn[2];

  char        *filename = "P2_triangle.vtk";
  int          n_vtx    = 7;
  PDM_g_num_t *vtx_g_num = malloc(sizeof(PDM_g_num_t) * 7);
  char        *vtx_field_name = "ho_bezier_basis";
  double      *vtx_field = malloc(sizeof(double) * 7);

  for (int i = 0; i < 7; i++) {
    vtx_g_num[i] = i + 1;
    vtx_field[i] = 0;
  }

  vtx_field[6] = 1;

  PDM_vtk_write_point_cloud_with_field(filename,
                                       n_vtx,
                                       tria_coord,
                                       vtx_g_num,
                                       NULL,
                                       1,
                      (const char **) &vtx_field_name,
                    (const double **) &vtx_field);

  PDM_MPI_Finalize ();

}
