#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_dbbtree.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_distrib.h"
#include "pdm_iso_surface.h"
#include "pdm_multipart.h"

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
_read_args(int                    argc,
           char                 **argv,
           PDM_g_num_t           *n_vtx_seg,
           PDM_Mesh_nodal_elt_t  *elt_type,
           double                *level,
           int                   *n_part,
           int                   *part_method)
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
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *level = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



static
inline
double
_unit_sphere
(
 double x,
 double y,
 double z
)
{
  return x * x + y * y + z * z - 0.125;
}

static
inline
void
_unit_sphere_gradient
(
 double  x,
 double  y,
 double  z,
 double *df_dx,
 double *df_dy,
 double *df_dz
)
{
  *df_dx = 2*x;
  *df_dy = 2*y;
  *df_dz = 2*z;
}


static
void
_generate_volume_mesh
(
 PDM_MPI_Comm                  comm,
 const PDM_g_num_t             n_vtx_seg,
 const PDM_Mesh_nodal_elt_t    elt_type,
 const PDM_split_dual_t        part_method,
 const int                     n_part,
 const double                  length,
       PDM_dcube_nodal_t     **_dcube,
       PDM_dmesh_nodal_t     **_dmn,
       PDM_multipart_t       **_mpart
)
{

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         -1.,
                                                         -1.,
                                                         -1.,
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  /*
   * Partitionnement
   */
  int n_zone = 1;
  // int n_part_zones = {n_part};
  int *n_part_zones = (int *) malloc(sizeof(int) * n_zone);
  n_part_zones[0] = n_part;
  printf("n_part = %d\n", n_part);
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

  *_dcube = dcube;
  *_dmn   = dmn;
  *_mpart = mpart;

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

  PDM_g_num_t          n_vtx_seg = 10;
  double               length    = 1.;
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_TETRA4;
  double                level     = 1.e-2;
  int                   n_part    = 1;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif
  //  5 -> tetra
  //  6 -> pyramid
  //  7 -> prism
  //  8 -> hexa
  //  9 -> poly3d

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &elt_type,
             &level,
             &n_part,
     (int *) &part_method);

  assert(PDM_Mesh_nodal_elt_dim_get(elt_type) == 3);

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
   *  Generate volume mesh
   */
  if (i_rank == 0) {
    printf("-- Generate volume mesh\n");
    fflush(stdout);
  }

  PDM_dcube_nodal_t     *dcube = NULL;
  PDM_dmesh_nodal_t     *dmn    = NULL;
  PDM_multipart_t       *mpart  = NULL;
  _generate_volume_mesh (comm,
                         n_vtx_seg,
                         elt_type,
                         part_method,
                         n_part,
                         length,
                         &dcube,
                         &dmn,
                         &mpart);

  // PDM_multipart_t* mpart_surf = _generate_surf_mesh (comm,
  //                                                    n_vtx_seg,
  //                                                    elt_type,
  //                                                    part_method,
  //                                                    n_part,
  //                                                    length);


  int          *pn_vtx         = (int *)          malloc(sizeof(int)           * n_part);
  int          *pn_cell        = (int *)          malloc(sizeof(int)           * n_part);
  double      **pvtx_coord     = (double **)      malloc(sizeof(double *)      * n_part);
  PDM_g_num_t **pvtx_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **pcell_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    pn_vtx[i_part] = PDM_multipart_part_vtx_coord_get(mpart,
                                                      0,
                                                      i_part,
                                                      &pvtx_coord[i_part],
                                                      PDM_OWNERSHIP_KEEP);

    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_VERTEX,
                                    &pvtx_ln_to_gn[i_part],
                                    PDM_OWNERSHIP_KEEP);

    pn_cell[i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                      0,
                                                      i_part,
                                                      PDM_MESH_ENTITY_CELL,
                                                      &pcell_ln_to_gn[i_part],
                                                      PDM_OWNERSHIP_KEEP);
  }


  PDM_part_mesh_nodal_elmts_t* pmne_vol = PDM_dmesh_nodal_to_part_mesh_nodal_elmts(dmn,
                                                                                   PDM_GEOMETRY_KIND_VOLUMIC,
                                                                                   n_part, // n_part
                                                                                   pn_vtx,
                                                                                   pvtx_ln_to_gn,
                                                                                   pn_cell,
                                                                                   pcell_ln_to_gn,
                                                                                   NULL);


  for (int i_part = 0; i_part < n_part; i_part++) {

    int id_section = 0;
    PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get(pmne_vol, id_section);
    int         *elmt_vtx                 = NULL;
    int         *parent_num               = NULL;
    PDM_g_num_t *numabs                   = NULL;
    PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
    PDM_part_mesh_nodal_elmts_block_std_get(pmne_vol, id_section, i_part, &elmt_vtx, &numabs, &parent_num, &parent_entitity_ln_to_gn);

    char filename[999];
    sprintf(filename, "out_volumic_%i_%i.vtk", i_part, i_rank);
    PDM_vtk_write_std_elements(filename,
                               pn_vtx[i_part],
                               pvtx_coord[i_part],
                               pvtx_ln_to_gn[i_part],
                               t_elt,
                               pn_cell[i_part],
                               elmt_vtx,
                               pcell_ln_to_gn[i_part],
                               0,
                               NULL,
                               NULL);
  }


  //

  int           n_part_mesh = 1;
  int          *part_n_elt       = malloc (sizeof(int          ) * n_part_mesh);
  double      **part_elt_extents = malloc (sizeof(double      *) * n_part_mesh);
  PDM_g_num_t **part_elt_g_num   = malloc (sizeof(PDM_g_num_t *) * n_part_mesh);

  double tolerance = 1.e-11;
  for(int i_part = 0; i_part < n_part_mesh; ++i_part) {
    part_elt_extents[i_part] = (double *) malloc( 6 * part_n_elt[i_part] * sizeof(double));

    double      *vtx_coord    = NULL;
    int         *face_vtx_idx = NULL;
    int         *face_vtx     = NULL;
    PDM_g_num_t *face_g_num   = NULL;

    double *_extents = part_elt_extents[i_part];
    for (int j = 0; j < part_n_elt[i_part]; j++) {

      for (int k1 = 0; k1 < 3; k1++) {
        _extents[  k1] = DBL_MAX;
        _extents[3+k1] = -DBL_MAX;
      }

      for (int k = face_vtx_idx[j]; k < face_vtx_idx[j+1]; k++) {
        int i_vtx = face_vtx[k] - 1;
        double *_coords = (double *) vtx_coord + 3 * i_vtx;

        for (int k1 = 0; k1 < 3; k1++) {
          _extents[k1]   = PDM_MIN (_coords[k1], _extents[k1]);
          _extents[3+k1] = PDM_MAX (_coords[k1], _extents[3+k1]);
        }
      }

      double delta = -DBL_MAX;

      for (int k1 = 0; k1 < 3; k1++) {
        delta = PDM_MAX (delta, fabs (_extents[k1+3] - _extents[k1]));
      }

      delta *= tolerance;

      for (int k1 = 0; k1 < 3; k1++) {
        _extents[k1]   +=  - delta;
        _extents[3+k1] +=    delta;
      }
      _extents += 6;
    }
  }




  /* Compute local extents */
  double my_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    for (int i = 0; i < part_n_elt[ipart]; i++) {
      for (int j = 0; j < 3; j++) {
        my_extents[j]   = PDM_MIN (my_extents[j],   part_elt_extents[ipart][6*i + j]);
        my_extents[j+3] = PDM_MAX (my_extents[j+3], part_elt_extents[ipart][6*i + 3 + j]);
      }
    }
  }



  double global_extents[6];
  PDM_MPI_Allreduce (my_extents,   global_extents,   3,
                     PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce (my_extents+3, global_extents+3, 3,
                     PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  /* Break symmetry */
  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX (max_range, global_extents[i+3] - global_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    global_extents[i]   -= max_range * 1.1e-3;
    global_extents[i+3] += max_range * 1.0e-3;
  }


  PDM_dbbtree_t *dbbt = PDM_dbbtree_create (comm, 3, global_extents);


  PDM_box_set_t  *surf_mesh_boxes = PDM_dbbtree_boxes_set (dbbt,
                                                           n_part_mesh,
                                                           part_n_elt,
                                      (const double **)    part_elt_extents,
                                  (const PDM_g_num_t **)   part_elt_g_num);


  int          n_line = 0;
  PDM_g_num_t *line_g_num = malloc(6 * n_line * sizeof(PDM_g_num_t));
  double      *line_coord = malloc(6 * n_line * sizeof(double     ));

  /*
   *  Construction des lignes
   */
  int          *pn_selected_vtx       = malloc(n_part * sizeof(int        ));
  PDM_g_num_t **pn_selected_vtx_g_num = malloc(n_part * sizeof(PDM_g_num_t));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    // int id_section = 0;
    // PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get(pmne_vol, id_section);
    // int         *elmt_vtx                 = NULL;
    // int         *parent_num               = NULL;
    // PDM_g_num_t *numabs                   = NULL;
    // PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
    // PDM_part_mesh_nodal_elmts_block_std_get(pmne_vol, id_section, i_part, &elmt_vtx, &numabs, &parent_num, &parent_entitity_ln_to_gn);

    pn_selected_vtx      [i_part] = 0;
    pn_selected_vtx_g_num[i_part] = malloc(pn_vtx[i_part] * sizeof(double));

    for(int i_vtx = 0; i_vtx < pn_vtx[i_part]; ++i_vtx) {

      int is_inside = 1;
      for (int j = 0; j < 3; j++) {
        if (pvtx_coord[i_part][3*i_vtx+j] < global_extents[j] ||
            pvtx_coord[i_part][3*i_vtx+j] > global_extents[j+3]) {
          is_inside = 0;
        }
      }

      if(is_inside == 1) {
        pn_selected_vtx_g_num[i_part][pn_selected_vtx[i_part]++] = pvtx_ln_to_gn[i_part][i_vtx];
      }

    }

    pn_selected_vtx_g_num[i_part] = realloc(pn_selected_vtx_g_num[i_part], pn_selected_vtx[i_part] * sizeof(PDM_g_num_t));

  }

  if(1 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_long(pn_selected_vtx_g_num[i_part], pn_selected_vtx[i_part], "pn_selected_vtx_g_num ::") ;
    }
  }


  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                                      pn_selected_vtx_g_num,
                                                      NULL,
                                                      pn_selected_vtx,
                                                      n_part,
                                                      comm);

  PDM_g_num_t *distrib_selected_vtx_idx = PDM_part_to_block_distrib_index_get(ptb);
  int dn_vtx_selected = PDM_part_to_block_n_elt_block_get (ptb);
  // PDM_g_num_t *block_g_num_a = PDM_part_to_block_block_gnum_get (ptb);

  double *blk_vtx_coord = NULL;
  PDM_part_to_block_exch (ptb,
                          3 * sizeof(double),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) pvtx_coord,
                          NULL,
                (void **) &blk_vtx_coord);


  /*
   *  On a selectinné pas mal de vertex maintenant on :
   *     -
   */
  // for(int i_vtx = 0; i_vtx < dn_vtx_selected; ++i_vtx) {

  //   global_extents
  //   blk_vtx_coord

  // }



  int         *box_idx   = NULL;
  PDM_g_num_t *box_g_num = NULL;
  PDM_dbbtree_lines_intersect_boxes(dbbt,
                                    n_line,
                                    line_g_num,
                                    line_coord,
                                    &box_idx,
                                    &box_g_num);

  free(blk_vtx_coord);
  PDM_part_to_block_free(ptb);

  // if (1) {
  //   for (int i = 0; i < n_line; i++) {
  //     log_trace("line "PDM_FMT_G_NUM": ", line_ln_to_gn[i]);
  //     for (int j = intersecting_box_idx[i]; j < intersecting_box_idx[i+1]; j++) {
  //       log_trace(PDM_FMT_G_NUM" ", intersecting_box_g_num[j]);
  //     }
  //     log_trace("\n");
  //   }
  // }



  free(line_g_num);
  free(line_coord);


  PDM_dbbtree_free (dbbt);
  PDM_box_set_destroy (&surf_mesh_boxes);



  for(int i_part = 0; i_part < n_part_mesh; ++i_part) {
    free(part_elt_extents);
    free(part_elt_g_num);
  }
  free(part_n_elt);

  PDM_part_mesh_nodal_elmts_free(pmne_vol);

  PDM_dcube_nodal_gen_free(dcube);
  PDM_multipart_free(mpart);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

}
