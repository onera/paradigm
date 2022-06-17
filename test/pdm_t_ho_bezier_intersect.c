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

  // Do triangle line intersection

  PDM_triangle_line_intersection(line, tria_coord, ip);

  // Get xyz coordinates from uv coordinates


  double *uvw     = malloc(sizeof(double) * 3 * 9);
  double *weights = malloc(sizeof(double) * 3 * 9);

  PDM_ho_bezier_basis(PDM_MESH_NODAL_TRIA3,
                      1,
                      3,
                      uvw,
                      weigths);


  PDM_ho_bezier_basis(PDM_MESH_NODAL_TRIAHO,
                      2,
                      6,
                      uvw,
                      weigths);

}
