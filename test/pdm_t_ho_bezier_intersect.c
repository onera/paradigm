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
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_ho_bezier_basis.h"
#include "pdm_triangle.h"

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

  // Set-up test
  double *line       = malloc(sizeof(double) * 6);
  double *tria_coord = malloc(sizeof(double) * 9);
  double *ip         = malloc(sizeof(double) * 3);

  line[0] = 0; line[3] = 0;
  line[1] = 0; line[4] = 0;
  line[2] = 0; line[5] = 0;

  tria_coord[0] = 0; tria_coord[3] = 0; tria_coord[6] = 0;
  tria_coord[1] = 0; tria_coord[4] = 0; tria_coord[7] = 0;
  tria_coord[2] = 0; tria_coord[5] = 0; tria_coord[8] = 0;

  // Do triangle P1 line intersection

  PDM_triangle_line_intersection(line, tria_coord, ip);

  // Construct matrices.  Since we have over determined system, need to find
  // which 2 out of 3 equations to use to develop equations. (Any 2 should
  // work since we've projected point to plane.)

  double rhs[2], c1[2], c2[2];

  for (int i = 0; i < 2; i++) {
    rhs[i] = ip[i] - tria_coord[6 + i];
    c1[i] = tria_coord[i] - tria_coord[6 + i];
    c2[i] = tria_coord[3 + i] - tria_coord[6 + i];
  }

  double det = PDM_DETERMINANT2X2(c1,c2);

  double pcoords[3];

  pcoords[0] = PDM_DETERMINANT2X2(rhs,c2) / det;
  pcoords[1] = PDM_DETERMINANT2X2(c1,rhs) / det;

  double weights[3];
  double closest_point[3];

  weights[0] = 1 - (pcoords[0] + pcoords[1]);
  weights[1] = pcoords[0];
  weights[2] = pcoords[1];
  if ( weights[0] >= 0.0 && weights[0] <= 1.0 &&
       weights[1] >= 0.0 && weights[1] <= 1.0 &&
       weights[2] >= 0.0 && weights[2] <= 1.0 ) {

    // Projection distance

    closest_point[0] = ip[0];
    closest_point[1] = ip[1];
    closest_point[2] = ip[2];

  } // end is inside case

  // Get uvw in triangle Pn

  double uv_clossest_Pn[3];

  double _uvPn_sub_tria[6];

  _uvPn_sub_tria[0] = uv_nodes[2*idx1];
  _uvPn_sub_tria[1] = uv_nodes[2*idx1+1];
  _uvPn_sub_tria[2] = uv_nodes[2*idx2];
  _uvPn_sub_tria[3] = uv_nodes[2*idx2+1];
  _uvPn_sub_tria[4] = uv_nodes[2*idx3];
  _uvPn_sub_tria[5] = uv_nodes[2*idx3+1];

  for (int j = 0; j < 2; j++) {
    for (int k = 0; k < 3; k++) {
      uv_clossest_Pn[j] += weights[k] * uvPn_sub_tria[2*k + j];
    }
  }

  // Get xyz coordinates from uvw coordinates

  double *uvw     = malloc(sizeof(double) * 3 * 9);
  double *weightsPn = malloc(sizeof(double) * 3 * 9);

  PDM_ho_bezier_basis(PDM_MESH_NODAL_TRIAHO,
                      2,
                      6,
                      uv_clossest_Pn,
                      weigthsPn);

}
