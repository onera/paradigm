
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_vtk.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Maximum number of sections depending of section type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/


void
PDM_vtk_write_boxes
(
 const char        *filename,
 const int          n_box,
 const double      *box_extents,
 const PDM_g_num_t *box_g_num
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nboxes\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*n_box);
  for (int i = 0; i < n_box; i++) {
    const double *e = box_extents + 6*i;
    fprintf(f, "%f %f %f\n", e[0], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[5]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[5]);
  }

  fprintf(f, "CELLS %d %d\n", n_box, 9*n_box);
  for (int i = 0; i < n_box; i++) {
    int j = 8*i;
    fprintf(f, "8 %d %d %d %d %d %d %d %d\n", j, j+1, j+2, j+3, j+4, j+5, j+6, j+7);
  }

  fprintf(f, "CELL_TYPES %d\n", n_box);
  for (int i = 0; i < n_box; i++) {
    fprintf(f, "12\n");
  }

  fprintf(f, "CELL_DATA %d\n", n_box);
  fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
  for (int i = 0; i < n_box; i++) {
    fprintf(f, ""PDM_FMT_G_NUM"\n", box_g_num[i]);
  }

  fclose(f);
};

void
PDM_vtk_write_circles
(
 const char        *filename,
 const int          n_circles,
 const double      *center,
 const double      *radius,
 const PDM_g_num_t *g_num,
 const int         *color,
 const int          resolution
)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "circles\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  double step = 2. * PDM_PI / (double) (resolution - 1);

  fprintf(f, "POINTS %d double\n", n_circles * resolution);
  for (int i = 0; i < n_circles; i++) {
    for (int j = 0; j < resolution; j++) {
      fprintf(f, "%f %f %f\n",
              center[3*i]     + radius[i]*cos(j*step),
              center[3*i + 1] + radius[i]*sin(j*step),
              center[3*i + 2]);
    }
  }

  fprintf(f, "CELLS %d %d\n", n_circles, n_circles * (resolution + 1));
  for (int i = 0; i < n_circles; i++) {
    fprintf(f, "%d \n", resolution);
    for (int j = 0; j < resolution-1; j++) {
      fprintf(f, "%d ", resolution*i + j);
    }
    fprintf(f, "%d\n", resolution*i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_circles);
  for (int i = 0; i < n_circles; i++) {
    fprintf(f, "4\n");
  }

  if (g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_circles);
    fprintf(f, "SCALARS gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_circles; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", g_num[i]);
    }
  }

  else if (color != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_circles);
    fprintf(f, "SCALARS color long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_circles; i++) {
      fprintf(f, "%d\n", color[i]);
    }
  }

  fclose(f);
}
