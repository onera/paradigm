#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_inside_cloud_surf.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dbbtree.h"
#include "pdm_surf_mesh.h"

#include "pdm_vtk.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

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
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

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
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_vtx_seg,
           PDM_g_num_t   *n_vtx_seg_cloud,
           PDM_g_num_t   *n_vtx_seg_tgt,
           double        *length,
           int           *post)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg_cloud = atol(argv[i]);
        *n_vtx_seg_cloud = (PDM_g_num_t) _n_vtx_seg_cloud;
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg_tgt = atol(argv[i]);
        *n_vtx_seg_tgt = (PDM_g_num_t) _n_vtx_seg_tgt;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static
void
_generate_cartesian_cloud
(
  int           n_vtx_x,
  int           n_vtx_y,
  int           n_vtx_z,
  double        length,
  double      **pt_coord_out,
  PDM_g_num_t **pt_gnum_out
)
{
  PDM_g_num_t n_pts = n_vtx_x * n_vtx_y * n_vtx_z;
  double      *pt_coord = malloc( 3 * n_pts * sizeof(double));
  PDM_g_num_t *pt_gnum  = malloc(     n_pts * sizeof(PDM_g_num_t));

  double step_x = length / (double) (n_vtx_x - 1);
  double step_y = length / (double) (n_vtx_y - 1);
  double step_z = length / (double) (n_vtx_z - 1);

  for (int i_vtx = 0; i_vtx < n_pts; ++i_vtx) {


    pt_gnum[i_vtx] = i_vtx;

    PDM_g_num_t indi = i_vtx % n_vtx_x;
    PDM_g_num_t indj = ((i_vtx - indi) / n_vtx_x) % n_vtx_y;
    PDM_g_num_t indk = i_vtx / (n_vtx_x * n_vtx_y);

    pt_coord[3 * i_vtx    ] = indi * step_x; //+ dcube->zero_x;
    pt_coord[3 * i_vtx + 1] = indj * step_y; //+ dcube->zero_y;
    pt_coord[3 * i_vtx + 2] = indk * step_z; //+ dcube->zero_z;

  }

  *pt_coord_out = pt_coord;
  *pt_gnum_out  = pt_gnum;

}


static
void
_generate_cartesian_boxes
(
  int           n_vtx_x,
  int           n_vtx_y,
  int           n_vtx_z,
  double        length,
  double      **box_coord_out,
  PDM_g_num_t **box_gnum_out
)
{
  PDM_g_num_t n_box = n_vtx_x * n_vtx_y * n_vtx_z;
  double      *box_coord = malloc( 6 * n_box * sizeof(double));
  PDM_g_num_t *box_gnum  = malloc(     n_box * sizeof(PDM_g_num_t));

  double step_x = length / (double) (n_vtx_x - 1);
  double step_y = length / (double) (n_vtx_y - 1);
  double step_z = length / (double) (n_vtx_z - 1);

  for (int i_box = 0; i_box < n_box; ++i_box) {


    box_gnum[i_box] = i_box;

    PDM_g_num_t indi = i_box % n_vtx_x;
    PDM_g_num_t indj = ((i_box - indi) / n_vtx_x) % n_vtx_y;
    PDM_g_num_t indk = i_box / (n_vtx_x * n_vtx_y);

    box_coord[6 * i_box    ] = indi * step_x; //+ dcube->zero_x;
    box_coord[6 * i_box + 1] = indj * step_y; //+ dcube->zero_y;
    box_coord[6 * i_box + 2] = indk * step_z; //+ dcube->zero_z;

    // PDM_g_num_t indi = i_box % n_vtx_x;
    // PDM_g_num_t indj = ((i_box - indi) / n_vtx_x) % n_vtx_y;
    // PDM_g_num_t indk = i_box / (n_vtx_x * n_vtx_y);

    box_coord[6 * i_box + 3] = (indi+1) * step_x; //+ dcube->zero_x;
    box_coord[6 * i_box + 4] = (indj+1) * step_y; //+ dcube->zero_y;
    box_coord[6 * i_box + 5] = (indk+1) * step_z; //+ dcube->zero_z;

  }

  *box_coord_out = box_coord;
  *box_gnum_out  = box_gnum;

}

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t           n_vtx_seg       = 10;
  PDM_g_num_t           n_vtx_seg_cloud = 10;
  PDM_g_num_t           n_vtx_seg_tgt   = 10;
  double                length          = 1.;
  int                   post            = 1;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &n_vtx_seg_cloud,
             &n_vtx_seg_tgt,
             &length,
             &post);

  double radius         = length;//2*length;

  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   *  Generate cloud
   */
  double      *pt_coord = NULL;
  PDM_g_num_t *pt_gnum  = NULL;
  PDM_g_num_t n_pts = n_vtx_seg_cloud * n_vtx_seg_cloud * n_vtx_seg_cloud;
  _generate_cartesian_cloud(n_vtx_seg_cloud,
                            n_vtx_seg_cloud,
                            n_vtx_seg_cloud,
                            length,
                            &pt_coord,
                            &pt_gnum);

  if(post) {
    PDM_vtk_write_point_cloud("pts_cloud.vtk", n_pts, pt_coord, pt_gnum, NULL);
  }

  /*
   *  Generate target_boxes
   */
  double      *tgt_box_extents = NULL;
  PDM_g_num_t *tgt_box_gnum    = NULL;
  PDM_g_num_t n_gtg_box = (n_vtx_seg_tgt - 1) * (n_vtx_seg_tgt - 1) * (n_vtx_seg_tgt - 1);
  _generate_cartesian_boxes(n_vtx_seg_tgt-1,
                            n_vtx_seg_tgt-1,
                            n_vtx_seg_tgt-1,
                            length,
                            &tgt_box_extents,
                            &tgt_box_gnum);
  if(post) {
    PDM_vtk_write_boxes("target_boxes.vtk", n_gtg_box, tgt_box_extents, tgt_box_gnum);
  }

  /*
   *  Generate boxes
   */
  double      *box_extents = NULL;
  PDM_g_num_t *box_gnum    = NULL;
  PDM_g_num_t n_box = (n_vtx_seg - 1) * (n_vtx_seg - 1) * (n_vtx_seg - 1);
  _generate_cartesian_boxes(n_vtx_seg-1,
                            n_vtx_seg-1,
                            n_vtx_seg-1,
                            length,
                            &box_extents,
                            &box_gnum);
  if(post) {
    PDM_vtk_write_boxes("boxes.vtk", n_box, box_extents, box_gnum);
  }


  /*
   * Start box tree
   */






  free(box_extents);
  free(box_gnum);
  free(tgt_box_extents);
  free(tgt_box_gnum);
  free(pt_coord);
  free(pt_gnum);


  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);
  if (i_rank == 0) {
  PDM_printf ("-- End\n");
  }

  PDM_MPI_Finalize ();

  return 0;
}
