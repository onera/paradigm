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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dmesh_partitioning.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dcube_gen.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_gnum.h"
#include "pdm_part1_to_selected_part2.h"
#include "pdm_part_to_block.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_vtk.h"
#include "pdm_dmesh.h"
#include "pdm_logging.h"
#include "pdm_priv.h"

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
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
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
 * \param [inout]   part_method Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length,
           int           *n_part,
           int           *post,
           int           *part_method)
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
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
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

  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_part    = 1;
  int                post      = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
     (int *) &part_method);

  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t  *dface_cell = NULL;
  int          *dface_vtx_idx = NULL;
  PDM_g_num_t  *dface_vtx = NULL;
  double       *dvtx_coord = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group = NULL;
  int           dface_vtxL;
  int           dFaceGroupL;

  PDM_dcube_t* dcube = PDM_dcube_gen_init(comm,
                                          n_vtx_seg,
                                          length,
                                          0.,
                                          0.,
                                          0.,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                        &n_face_group,
                        &dn_cell,
                        &dn_face,
                        &dn_vtx,
                        &dface_vtxL,
                        &dFaceGroupL);

  PDM_dcube_gen_data_get(dcube,
                         &dface_cell,
                         &dface_vtx_idx,
                         &dface_vtx,
                         &dvtx_coord,
                         &dface_group_idx,
                         &dface_group);

  /*
   * Create dmesh
   */
  int n_jn = 0;
  PDM_dmesh_t* dm = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                     dn_cell,
                                     dn_face,
                                     -1, // dn_edge
                                     dn_vtx,
                                     n_face_group, // n_bnd
                                     n_jn,
                                     comm);

  int *dface_join_idx = (int *) malloc( (n_jn+1) * sizeof(int));
  dface_join_idx[0] = 0;
  PDM_dmesh_set(dm,
                dvtx_coord,
                dface_vtx_idx,
                dface_vtx,
                dface_cell,
                dface_group_idx,
                dface_group,
                NULL,
                dface_join_idx,
                NULL);

  /*
   * Partitionnement
   */
  int n_zone = 1;
  int n_part_zones = n_part;
  int mpart_id = PDM_multipart_create(n_zone,
                                      &n_part_zones,
                                      PDM_FALSE,
                                      part_method,
                                      PDM_PART_SIZE_HOMOGENEOUS,
                                      NULL,
                                      comm,
                                      PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart_id, -1, "PDM_PART_RENUM_CELL_NONE", NULL, "PDM_PART_RENUM_FACE_NONE");
  PDM_multipart_register_block(mpart_id, 0, dm);
  PDM_multipart_run_ppart(mpart_id);

  /*
   * Get the partition zone
   */
  int i_zone = 0;

  double      **cell_center        = (double      **) malloc( n_part_zones * sizeof(double      *));
  PDM_g_num_t **selected_g_num     = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pcell_ln_to_gn     = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pface_ln_to_gn     = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  int         **selected_g_num_idx = (int         **) malloc( n_part_zones * sizeof(int         *));
  int          *pn_cell            = (int          *) malloc( n_part_zones * sizeof(int          ));
  int          *pn_select_cell     = (int          *) malloc( n_part_zones * sizeof(int          ));
  double      **weight             = (double      **) malloc( n_part_zones * sizeof(double      *));
  int         **pcell_face         = (int         **) malloc( n_part_zones * sizeof(int         *));
  int         **pcell_face_idx     = (int         **) malloc( n_part_zones * sizeof(int         *));

  PDM_gen_gnum_t* gnum_extract = PDM_gnum_create(3,
                                                 n_part_zones,
                                                 PDM_FALSE,
                                                 1.e-6,
                                                 comm,
                                                 PDM_OWNERSHIP_KEEP);

  for (int i_part = 0; i_part < n_part_zones; i_part++){

    int n_proc, tn_part;
    int n_cell, n_face, n_vtx, n_bounds, n_joins, n_part_joins;
    int scell_face, sface_vtx, sface_bound, sface_join;
    int  n_section;
    int* n_elt;

    PDM_multipart_part_dim_get(mpart_id, i_zone, i_part, &n_section, &n_elt,
                               &n_cell, &n_face, &n_part_joins, &n_vtx, &n_proc, &tn_part,
                               &scell_face, &sface_vtx, &sface_bound, &n_bounds, &sface_join, &n_joins);

    double       *vtx;
    int          *cell_face_idx, *cell_face, *face_cell, *face_vtx_idx, *face_vtx;
    int          *face_bound_idx, *face_bound, *face_join_idx, *face_join;
    int          *face_part_bound_proc_idx, *face_part_bound_part_idx, *face_part_bound;
    PDM_g_num_t  *cell_ln_to_gn, *face_ln_to_gn, *vtx_ln_to_gn, *face_bound_ln_to_gn, *face_join_ln_to_gn;
    int          *cell_tag, *face_tag, *vtx_tag;
    int         **elt_vtx_idx;
    int         **elt_vtx;
    PDM_g_num_t **elt_section_ln_to_gn;


    PDM_multipart_part_val_get(mpart_id, i_zone, i_part, &elt_vtx_idx, &elt_vtx, &elt_section_ln_to_gn,
                               &cell_tag, &cell_face_idx, &cell_face, &cell_ln_to_gn,
                               &face_tag, &face_cell, &face_vtx_idx, &face_vtx, &face_ln_to_gn,
                               &face_part_bound_proc_idx, &face_part_bound_part_idx, &face_part_bound,
                               &vtx_tag, &vtx, &vtx_ln_to_gn, &face_bound_idx, &face_bound,
                               &face_bound_ln_to_gn, &face_join_idx, &face_join, &face_join_ln_to_gn);

    pn_cell       [i_part] = n_cell;
    pcell_ln_to_gn[i_part] = cell_ln_to_gn;
    pface_ln_to_gn[i_part] = face_ln_to_gn;
    pcell_face    [i_part] = cell_face;
    pcell_face_idx[i_part] = cell_face_idx;


    /*
     * Compute center-cell and extract cells corresponding to criteria
     */
    double *face_center         = (double *) malloc( 3 * n_face * sizeof(double));

    for(int i_face = 0; i_face < n_face; ++i_face) {
      face_center[3*i_face  ] = 0.;
      face_center[3*i_face+1] = 0.;
      face_center[3*i_face+2] = 0.;
      int n_vtx_on_face = face_vtx_idx[i_face+1] - face_vtx_idx[i_face];
      for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = PDM_ABS(face_vtx[idx_vtx])-1;
        face_center[3*i_face  ] += vtx[3*i_vtx  ];
        face_center[3*i_face+1] += vtx[3*i_vtx+1];
        face_center[3*i_face+2] += vtx[3*i_vtx+2];
      }
      face_center[3*i_face  ] = face_center[3*i_face  ] / n_vtx_on_face;
      face_center[3*i_face+1] = face_center[3*i_face+1] / n_vtx_on_face;
      face_center[3*i_face+2] = face_center[3*i_face+2] / n_vtx_on_face;
    }

    cell_center[i_part] = (double *) malloc( 3 * n_cell * sizeof(double));

    for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

      cell_center[i_part][3*i_cell  ] = 0.;
      cell_center[i_part][3*i_cell+1] = 0.;
      cell_center[i_part][3*i_cell+2] = 0.;

      int n_face_on_cell = cell_face_idx[i_cell+1] - cell_face_idx[i_cell];

      for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {

        int i_face = PDM_ABS(cell_face[idx_face])-1;
        cell_center[i_part][3*i_cell  ] += face_center[3*i_face  ];
        cell_center[i_part][3*i_cell+1] += face_center[3*i_face+1];
        cell_center[i_part][3*i_cell+2] += face_center[3*i_face+2];
      }
      cell_center[i_part][3*i_cell  ] = cell_center[i_part][3*i_cell  ] / n_face_on_cell;
      cell_center[i_part][3*i_cell+1] = cell_center[i_part][3*i_cell+1] / n_face_on_cell;
      cell_center[i_part][3*i_cell+2] = cell_center[i_part][3*i_cell+2] / n_face_on_cell;
    }

    free(face_center);

    selected_g_num    [i_part] = (PDM_g_num_t *) malloc(  n_cell      * sizeof(PDM_g_num_t));
    selected_g_num_idx[i_part] = (int         *) malloc( (n_cell + 1) * sizeof(int        ));

    /*
     * Sub-part
     */
    double bbox[6];
    bbox[0] = 0.2;
    bbox[1] = 0.2;
    bbox[2] = 0.2;
    bbox[3] = 0.8;
    bbox[4] = 0.8;
    bbox[5] = 0.8;
    int n_select_cell = 0;
    selected_g_num_idx[i_part][0] = 0;
    for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

      int inside = 1;
      for(int i = 0; i < 3; ++i) {
        if (cell_center[i_part][3*i_cell+i] > bbox[i+3] || cell_center[i_part][3*i_cell+i] < bbox[i]) {
          inside = 0;
        }
      }
      if(inside == 1) {
        selected_g_num[i_part][n_select_cell++] = cell_ln_to_gn[i_cell];
      }
      selected_g_num_idx[i_part][i_cell+1] = selected_g_num_idx[i_part][i_cell] + inside;
    }

    selected_g_num[i_part] = realloc(selected_g_num[i_part], n_select_cell * sizeof(PDM_g_num_t));
    pn_select_cell[i_part] = n_select_cell;

    weight        [i_part] = malloc(n_select_cell * sizeof(double));
    for(int i = 0; i < n_select_cell; ++i) {
      weight        [i_part][i] = 1.;
    }

    PDM_log_trace_array_long(selected_g_num    [i_part], n_select_cell, "selected_g_num     : ");
    PDM_log_trace_array_int (selected_g_num_idx[i_part], n_cell+1     , "selected_g_num_idx : ");
    PDM_log_trace_array_long(cell_ln_to_gn             , n_cell       , "cell_ln_to_gn      : ");


  }
  free(dface_join_idx);

  /*
   *  We know the extraction of cells required : now we use part_to_part to resetup a coherent part :
   *           - Redistribute over all process ( And order with scoth or something else)
   *           - Keep parent_cell_g_num
   *           - Rebuild by descending connectivity a coherent partition (cell -> face -> vtx + face_group)
   */

  PDM_g_num_t **child_selected_g_num = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  for (int i_part = 0; i_part < n_part_zones; i_part++){
    PDM_gnum_set_from_coords(gnum_extract, i_part, pn_cell[i_part], cell_center[i_part], NULL);
  }

  PDM_gnum_compute(gnum_extract);

  for (int i_part = 0; i_part < n_part_zones; i_part++){
    child_selected_g_num[i_part] = PDM_gnum_get(gnum_extract, i_part);
  }

  /*
   *  Remake equilibrate block -> Block is not partial
   */
  PDM_part_to_block_t *ptb_equi = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           child_selected_g_num,
                                                           weight,
                                                           pn_select_cell,
                                                           n_part_zones,
                                                           comm);

  int dn_cell_equi = PDM_part_to_block_n_elt_block_get (ptb_equi);
  PDM_g_num_t *dextract_gnum = PDM_part_to_block_block_gnum_get(ptb_equi);


  PDM_part1_to_selected_part2_t* ptp = PDM_part1_to_selected_part2_create((const PDM_g_num_t **) pcell_ln_to_gn,
                                                                          pn_cell,
                                                                          n_part_zones,
                                                                          (const PDM_g_num_t **) &dextract_gnum,
                                                                          &dn_cell_equi,
                                                                          1,
                                                                          (const int         **) selected_g_num_idx,
                                                                          (const PDM_g_num_t **) child_selected_g_num,
                                                                          comm);

  /*
   * Extract cell_face
   */
  int         **pextract_cell_face_idx = (int         **) malloc( n_part_zones * sizeof(int         *));
  PDM_g_num_t **pextract_cell_face     = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  for (int i_part = 0; i_part < n_part_zones; i_part++){

    pextract_cell_face_idx[i_part] = malloc( (pn_select_cell[i_part]+1) * sizeof(int));

    int i_extract_cell = 0;
    pextract_cell_face_idx[i_part][0] = 0;
    for(int i_cell = 0; i_cell < pn_cell[i_part]; ++i_cell) {
      for(int s_cell = selected_g_num_idx[i_part][i_cell]; s_cell < selected_g_num_idx[i_part][i_cell+1]; ++s_cell) {
        pextract_cell_face_idx[i_part][i_extract_cell+1] = pextract_cell_face_idx[i_part][i_extract_cell] + (pcell_face_idx[i_part][i_cell+1] - pcell_face_idx[i_part][i_cell]);
        i_extract_cell++;
      }
    }

    assert(i_extract_cell == pn_select_cell[i_part]);

    pextract_cell_face[i_part] = malloc( (pextract_cell_face_idx[i_part][pn_select_cell[i_part]]) * sizeof(PDM_g_num_t));

    int idx_write = 0;
    for(int i_cell = 0; i_cell < pn_cell[i_part]; ++i_cell) {
      for(int s_cell = selected_g_num_idx[i_part][i_cell]; s_cell < selected_g_num_idx[i_part][i_cell+1]; ++s_cell) {
        for(int idx_face = pcell_face_idx[i_part][i_cell]; idx_face < pcell_face_idx[i_part][i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS (pcell_face[i_part][idx_face])-1;
          int sgn    = PDM_SIGN(pcell_face[i_part][idx_face]);
          PDM_g_num_t g_face = pface_ln_to_gn[i_part][i_face];
          pextract_cell_face[i_part][idx_write++] = sgn * g_face;
        }
      }
    }


    PDM_log_trace_array_long(pextract_cell_face    [i_part], pextract_cell_face_idx[i_part][pn_select_cell[i_part]], "pextract_cell_face     : ");
    PDM_log_trace_array_int (pextract_cell_face_idx[i_part], pn_select_cell[i_part]+1                              , "pextract_cell_face_idx : ");


  }



  /*
   *  Exchange cell_face in global numebering = cell_face + face_ln_to_gn
   */





  PDM_part1_to_selected_part2_free(ptp);

  PDM_part_to_block_free(ptb_equi);
  PDM_gnum_free(gnum_extract);


  for (int i_part = 0; i_part < n_part_zones; i_part++){
    free(cell_center       [i_part]);
    free(selected_g_num    [i_part]);
    free(selected_g_num_idx[i_part]);
    free(weight[i_part]);
    free(pextract_cell_face[i_part]);
    free(pextract_cell_face_idx[i_part]);
  }
  free(cell_center);
  free(selected_g_num);
  free(selected_g_num_idx);
  free(pn_cell);
  free(pcell_face);
  free(pcell_face_idx);
  free(child_selected_g_num);
  free(pn_select_cell);
  free(weight);
  free(pcell_ln_to_gn);
  free(pface_ln_to_gn);
  free(pextract_cell_face);
  free(pextract_cell_face_idx);

  PDM_multipart_free(mpart_id);
  PDM_dcube_gen_free(dcube);
  PDM_dmesh_free(dm);

  PDM_MPI_Finalize();
  return 0;
}
