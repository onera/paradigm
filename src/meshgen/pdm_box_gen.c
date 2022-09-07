/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <sys/resource.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_timer.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"

#include "pdm_box_gen.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/


/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/* MPI ??? */
void
PDM_box_gen_cartesian
(
  int           n_vtx_x,
  int           n_vtx_y,
  int           n_vtx_z,
  double        length,
  int          *n_box_out,
  double      **box_coord_out,
  PDM_g_num_t **box_gnum_out
)
{
  PDM_g_num_t n_box = (n_vtx_x - 1) * ( n_vtx_y - 1) * ( n_vtx_z - 1);
  // PDM_g_num_t n_box = (n_vtx_x) * ( n_vtx_y) * ( n_vtx_z);
  *n_box_out = n_box;
  double      *box_coord = malloc( 6 * n_box * sizeof(double));
  PDM_g_num_t *box_gnum  = malloc(     n_box * sizeof(PDM_g_num_t));

  double step_x = length / (double) (n_vtx_x - 1);
  double step_y = length / (double) (n_vtx_y - 1);
  double step_z = length / (double) (n_vtx_z - 1);

  int n_box_x = n_vtx_x - 1;
  int n_box_y = n_vtx_y - 1;
  // int n_box_z = n_vtx_z - 1;

  for (int i_box = 0; i_box < n_box; ++i_box) {


    box_gnum[i_box] = i_box;

    PDM_g_num_t ind_box_i = i_box % n_box_x;
    PDM_g_num_t ind_box_j = ((i_box - ind_box_i) / n_box_x) % n_box_y;
    PDM_g_num_t ind_box_k = i_box / (n_box_x * n_box_y);

    int i_vtx = ind_box_i + ind_box_j * n_vtx_x + ind_box_k * n_vtx_x * n_vtx_y;

    PDM_g_num_t indi = i_vtx % n_vtx_x;
    PDM_g_num_t indj = ((i_vtx - indi) / n_vtx_x) % n_vtx_y;
    PDM_g_num_t indk = i_vtx / (n_vtx_x * n_vtx_y);

    box_coord[6 * i_box    ] = indi * step_x; //+ dcube->zero_x;
    box_coord[6 * i_box + 1] = indj * step_y; //+ dcube->zero_y;
    box_coord[6 * i_box + 2] = indk * step_z; //+ dcube->zero_z;

    box_coord[6 * i_box + 3] = (indi+1) * step_x; //+ dcube->zero_x;
    box_coord[6 * i_box + 4] = (indj+1) * step_y; //+ dcube->zero_y;
    box_coord[6 * i_box + 5] = (indk+1) * step_z; //+ dcube->zero_z;

  }

  *box_coord_out = box_coord;
  *box_gnum_out  = box_gnum;

}



/**
 *
 * \brief Generate a random set of boxes
 *
 * \param [in]   comm                   MPI Communicator id
 * \param [in]   seed                   Random seed
 * \param [in]   geometric_g_num        Compute global ids from coordinates
 * \param [in]   gn_box                 Global number of boxes
 * \param [in]   min_size               Minimal box size
 * \param [in]   max_size               Maximal box size
 * \param [in]   x_min                  Minimal X-coordinate for box centers
 * \param [in]   y_min                  Minimal Y-coordinate for box centers
 * \param [in]   z_min                  Minimal Z-coordinate for box centers
 * \param [in]   x_max                  Maximal X-coordinate for box centers
 * \param [in]   y_max                  Maximal Y-coordinate for box centers
 * \param [in]   z_max                  Maximal Z-coordinate for box centers
 * \param [out]  n_box                  Local number of boxes
 * \param [out]  box_extents            Extents of the local boxes
 * \param [out]  box_ln_to_gn           Global ids of the local boxes
 *
 */

void
PDM_box_gen_random
(
 PDM_MPI_Comm   comm,
 int            seed,
 int            geometric_g_num,
 PDM_g_num_t    gn_box,
 double         min_size,
 double         max_size,
 double         x_min,
 double         y_min,
 double         z_min,
 double         x_max,
 double         y_max,
 double         z_max,
 int           *n_box,
 double       **box_extents,
 PDM_g_num_t  **box_ln_to_gn
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  double origin[3] = {x_min, y_min, z_min};
  double length[3] = {x_max - x_min, y_max - y_min, z_max - z_min};

  /*
   *  Generate random boxes
   */
  PDM_g_num_t *distrib_box = PDM_compute_uniform_entity_distribution (comm,
                                                                      gn_box);
  *n_box = (int) (distrib_box[i_rank+1] - distrib_box[i_rank]);

  *box_extents = malloc (sizeof(double) * (*n_box) * 6);
  for (int i = 0; i < (*n_box); i++) {

    unsigned int _seed = (unsigned int) (distrib_box[i_rank] + i) + 1;
    srand(_seed + seed);

    for (int j = 0; j < 3; j++) {
      double mid = origin[j] + length[j] * ((double) rand() / (double) RAND_MAX);
      double size = min_size + 0.5*(max_size - min_size) * ((double) rand() / (double) RAND_MAX);

      (*box_extents)[6*i + j]     = mid - size;
      (*box_extents)[6*i + j + 3] = mid + size;
    }
  }


  if (geometric_g_num) {
    double *box_centers = malloc (sizeof(double) * (*n_box) * 3);
    for (int i = 0; i < (*n_box); i++) {
      for (int j = 0; j < 3; j++) {
        box_centers[3*i+j] = 0.5 * ((*box_extents)[6*i+j] + (*box_extents)[6*i+j+3]);
      }
    }
    PDM_gen_gnum_t *gen_gnum = PDM_gnum_create (3,
                                                1,
                                                PDM_FALSE,
                                                1.e-3,
                                                comm,
                                                PDM_OWNERSHIP_USER);

    PDM_gnum_set_from_coords (gen_gnum,
                              0,
                              *n_box,
                              box_centers,
                              NULL);

    PDM_gnum_compute (gen_gnum);
    free (box_centers);

    *box_ln_to_gn = PDM_gnum_get (gen_gnum, 0);

    PDM_gnum_free (gen_gnum);
  }
  else {
    *box_ln_to_gn = malloc(sizeof(PDM_g_num_t) * (*n_box));
    for (int i = 0; i < (*n_box); i++) {
      (*box_ln_to_gn)[i] = distrib_box[i_rank] + i + 1;
    }
  }
  free (distrib_box);

}
