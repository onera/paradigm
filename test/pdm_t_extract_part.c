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
#include "pdm_part_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_extract_part.h"
#include "pdm_vtk.h"
#include "pdm_dmesh.h"
#include "pdm_unique.h"
#include "pdm_part_geom.h"
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


static
void
PDM_pconnectivity_to_pconnectivity
(
  const PDM_MPI_Comm    comm,
  const int             n_part1,
  const int            *n_part1_entity1,
  const int           **part1_entity1_entity2_idx,
  const int           **part1_entity1_entity2,
  const PDM_g_num_t   **part1_entity1_ln_to_gn,
  const PDM_g_num_t   **part1_entity2_ln_to_gn,
  const int             n_part2,
  const int            *n_part2_entity1,
  const PDM_g_num_t   **part2_entity1_ln_to_gn,
  const int           **part2_entity1_to_part1_entity1_idx,
  const PDM_g_num_t   **part2_entity1_to_part1_entity1,
        int           **n_part2_entity2,
        int          ***part2_entity1_entity2_idx,
        int          ***part2_entity1_entity2,
        PDM_g_num_t  ***part2_entity2_ln_to_gn
)
{
  PDM_UNUSED(n_part2_entity2);
  PDM_UNUSED(part2_entity1_entity2_idx);
  PDM_UNUSED(part2_entity1_entity2);
  PDM_UNUSED(part2_entity2_ln_to_gn);

  PDM_part_to_part_t* ptp = PDM_part_to_part_create(part2_entity1_ln_to_gn,
                                                    n_part2_entity1,
                                                    n_part2,
                                                    part1_entity1_ln_to_gn,
                                                    n_part1_entity1,
                                                    n_part1,
                                                    part2_entity1_to_part1_entity1_idx,
                                                    part2_entity1_to_part1_entity1,
                                                    comm);

  /*
   * Protocol are created then we can extract information in part1 to reverse send it to part2
   */
  int          *n_ref_entity1     = NULL;
  int         **ref_l_num_entity1 = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_entity1, &ref_l_num_entity1);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);

  /* Create buffer */
  int         **send_entity1_entity2_n = malloc(n_part1 * sizeof(int         *));
  PDM_g_num_t **send_entity1_entity2   = malloc(n_part1 * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part1; ++i_part) {

    /*
     * Compute stride size
     */
    send_entity1_entity2_n[i_part] = malloc( gnum1_come_from_idx[i_part][n_ref_entity1[i_part]] * sizeof(int));

    int n_tot_send = 0;
    for(int j = 0; j < n_ref_entity1[i_part]; ++j) {
      for(int k = gnum1_come_from_idx[i_part][j]; k < gnum1_come_from_idx[i_part][j+1]; ++k) {
        int l_entity1     = ref_l_num_entity1[i_part][k]-1;
        int n_loc_entity2 = part1_entity1_entity2_idx[i_part][l_entity1+1] - part1_entity1_entity2_idx[i_part][l_entity1];
        send_entity1_entity2_n[i_part][k] = n_loc_entity2;
        n_tot_send += n_loc_entity2;
      }
    }

    // int* send_entity1_entity2_idx = malloc( (gnum1_come_from_idx[i_part][n_ref_entity1[i_part]] + 1) * sizeof(int));
    // send_entity1_entity2_idx[0] = 0;
    // for(int i = 0; i < gnum1_come_from_idx[i_part][n_ref_entity1[i_part]]; ++i) {
    //   send_entity1_entity2_idx[i+1] = send_entity1_entity2_idx[i] + send_entity1_entity2_n[i_part][i];
    // }

    send_entity1_entity2[i_part] = malloc( n_tot_send * sizeof(PDM_g_num_t));
    int idx_write = 0;
    for(int j = 0; j < n_ref_entity1[i_part]; ++j) {
      for(int k = gnum1_come_from_idx[i_part][j]; k < gnum1_come_from_idx[i_part][j+1]; ++k) {
        int l_face = ref_l_num_entity1[i_part][k]-1;
        for(int l = part1_entity1_entity2_idx[i_part][l_face]; l < part1_entity1_entity2_idx[i_part][l_face+1]; ++l) {
          send_entity1_entity2[i_part][idx_write++] = part1_entity2_ln_to_gn[i_part][part1_entity1_entity2[i_part][l]-1];
        }
      }
    }

    // printf("idx_write = %i | 4 * n_extract_face = %i \n", idx_write, 4 * n_extract_face);
    // PDM_log_trace_array_long(send_face_vtx[i_part], 4 * gnum1_come_from_idx[i_part][n_ref_face[i_part]], "send_face_vtx      : ");
  }

  int         **recv_entity1_entity2_n = NULL;
  PDM_g_num_t **recv_entity1_entity2   = NULL;
  int           exch_request = -1;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 -1,
                                 sizeof(PDM_g_num_t),
                (const int  **)  send_entity1_entity2_n,
                (const void **)  send_entity1_entity2,
                                 &recv_entity1_entity2_n,
                    (void ***)   &recv_entity1_entity2,
                                 &exch_request);

  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  for(int i_part = 0; i_part < n_part1; ++i_part) {
    free(send_entity1_entity2_n[i_part]);
    free(send_entity1_entity2  [i_part]);
  }
  free(send_entity1_entity2_n);
  free(send_entity1_entity2  );

  /*
   * Post-treatment
   */
  int          *_n_part2_entity2           = malloc(n_part2 * sizeof(int          ));
  int         **_part2_entity1_entity2_idx = malloc(n_part2 * sizeof(int         *));
  int         **_part2_entity1_entity2     = malloc(n_part2 * sizeof(int         *));
  PDM_g_num_t **_part2_entity2_ln_to_gn    = malloc(n_part2 * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part2; ++i_part) {

    _part2_entity1_entity2_idx[i_part] = malloc( (n_part2_entity1[i_part] + 1) * sizeof(int));

    /* Compute recv stride */
    _part2_entity1_entity2_idx[i_part][0] = 0;
    for(int i_entity1 = 0; i_entity1 < n_part2_entity1[i_part]; ++i_entity1) {
      _part2_entity1_entity2_idx[i_part][i_entity1+1] = _part2_entity1_entity2_idx[i_part][i_entity1] + recv_entity1_entity2_n[i_part][i_entity1];
    }
    int n_recv_entity1_entity2 = _part2_entity1_entity2_idx[i_part][n_part2_entity1[i_part]];

    _part2_entity2_ln_to_gn[i_part] = malloc( n_recv_entity1_entity2      * sizeof(PDM_g_num_t));

    int *unique_order_entity2     = (int         * ) malloc(n_recv_entity1_entity2 * sizeof(int        ));
    for(int i = 0; i < n_recv_entity1_entity2; ++i) {
      _part2_entity2_ln_to_gn[i_part][i] = PDM_ABS(recv_entity1_entity2[i_part][i]);
    }

    int n_extract_entity2 = PDM_inplace_unique_long2(_part2_entity2_ln_to_gn[i_part], unique_order_entity2, 0, n_recv_entity1_entity2-1);
    _part2_entity2_ln_to_gn[i_part] = realloc(_part2_entity2_ln_to_gn[i_part],  n_extract_entity2      * sizeof(PDM_g_num_t));

    /* Recompute local numbering */
    _part2_entity1_entity2 [i_part] = malloc( n_recv_entity1_entity2 * sizeof(int        ));

    for(int idx = 0; idx < n_recv_entity1_entity2; ++idx) {
      int g_sgn  = PDM_SIGN(recv_entity1_entity2[i_part][idx]);
      int l_elmt = unique_order_entity2[idx];
      _part2_entity1_entity2[i_part][idx] = (l_elmt + 1) * g_sgn;
    }
    free(unique_order_entity2);
  }

  for(int i_part = 0; i_part < n_part2; ++i_part) {
    free(recv_entity1_entity2_n[i_part]);
    free(recv_entity1_entity2  [i_part]);
  }
  free(recv_entity1_entity2_n);
  free(recv_entity1_entity2  );

  PDM_part_to_part_free(ptp);

  *n_part2_entity2           = _n_part2_entity2;
  *part2_entity1_entity2_idx = _part2_entity1_entity2_idx;
  *part2_entity1_entity2     = _part2_entity1_entity2;
  *part2_entity2_ln_to_gn    = _part2_entity2_ln_to_gn;

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
  PDM_multipart_t *mpart_id = PDM_multipart_create(n_zone,
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

  double      **cell_center             = (double      **) malloc( n_part_zones * sizeof(double      *));
  PDM_g_num_t **selected_g_num          = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  int         **selected_l_num          = (int         **) malloc( n_part_zones * sizeof(int         *));
  PDM_g_num_t **pcell_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pface_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pvtx_ln_to_gn           = (PDM_g_num_t **) malloc( n_part_zones * sizeof(PDM_g_num_t *));
  int         **selected_g_num_idx      = (int         **) malloc( n_part_zones * sizeof(int         *));
  int          *pn_cell                 = (int          *) malloc( n_part_zones * sizeof(int          ));
  int          *pn_face                 = (int          *) malloc( n_part_zones * sizeof(int          ));
  int          *pn_vtx                  = (int          *) malloc( n_part_zones * sizeof(int          ));
  int          *pn_select_cell          = (int          *) malloc( n_part_zones * sizeof(int          ));
  // double      **weight                  = (double      **) malloc( n_part_zones * sizeof(double      *));
  int         **pcell_face              = (int         **) malloc( n_part_zones * sizeof(int         *));
  int         **pcell_face_idx          = (int         **) malloc( n_part_zones * sizeof(int         *));
  int         **pface_vtx               = (int         **) malloc( n_part_zones * sizeof(int         *));
  int         **pface_vtx_idx           = (int         **) malloc( n_part_zones * sizeof(int         *));
  double      **pvtx_coord              = (double      **) malloc( n_part_zones * sizeof(double      *));
  // double      **tmp_extract_cell_center = (double      **) malloc( n_part_zones * sizeof(double      *));

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
    pvtx_ln_to_gn [i_part] = vtx_ln_to_gn;
    pcell_face    [i_part] = cell_face;
    pcell_face_idx[i_part] = cell_face_idx;
    pn_face       [i_part] = n_face;
    pn_vtx        [i_part] = n_vtx;

    pface_vtx    [i_part] = face_vtx;
    pface_vtx_idx[i_part] = face_vtx_idx;
    pvtx_coord   [i_part] = vtx;

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

    selected_g_num         [i_part] = (PDM_g_num_t *) malloc(  n_cell          * sizeof(PDM_g_num_t));
    selected_l_num         [i_part] = (int         *) malloc(  n_cell          * sizeof(int        ));
    selected_g_num_idx     [i_part] = (int         *) malloc( (n_cell + 1)     * sizeof(int        ));
    // tmp_extract_cell_center[i_part] = (double      *) malloc(  3 * n_cell      * sizeof(double     ));

    /*
     * Sub-part
     */

    double bbox[6];
    bbox[0] = 0.3;
    bbox[1] = 0.3;
    bbox[2] = 0.35;
    bbox[3] = 0.7;
    bbox[4] = 0.7;
    bbox[5] = 0.65;

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
        selected_g_num         [i_part][n_select_cell]     = cell_ln_to_gn[i_cell];
        selected_l_num         [i_part][n_select_cell]     = i_cell;
        // tmp_extract_cell_center[i_part][3*n_select_cell  ] = cell_center[i_part][3*i_cell  ];
        // tmp_extract_cell_center[i_part][3*n_select_cell+1] = cell_center[i_part][3*i_cell+1];
        // tmp_extract_cell_center[i_part][3*n_select_cell+2] = cell_center[i_part][3*i_cell+2];
        n_select_cell++;

      }
      selected_g_num_idx[i_part][i_cell+1] = selected_g_num_idx[i_part][i_cell] + inside;
    }

    selected_g_num         [i_part] = realloc(selected_g_num[i_part], n_select_cell * sizeof(PDM_g_num_t));
    selected_l_num         [i_part] = realloc(selected_l_num[i_part], n_select_cell * sizeof(int        ));
    pn_select_cell         [i_part] = n_select_cell;
    // tmp_extract_cell_center[i_part] = realloc(tmp_extract_cell_center[i_part], 3 * n_select_cell * sizeof(double));

    // PDM_log_trace_array_long(selected_g_num    [i_part], n_select_cell, "selected_g_num     : ");
    // PDM_log_trace_array_int (selected_g_num_idx[i_part], n_cell+1     , "selected_g_num_idx : ");
    // PDM_log_trace_array_long(cell_ln_to_gn             , n_cell       , "cell_ln_to_gn      : ");


  }
  free(dface_join_idx);

  /*
   * Extract
   */
  int n_part_out = 1;
  PDM_extract_part_t* extrp = PDM_extract_part_create(3,
                                                      n_part,
                                                      n_part_out,
                                                      PDM_SPLIT_DUAL_WITH_HILBERT,
                                                      PDM_OWNERSHIP_KEEP,
                                                      comm);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_extract_part_part_set(extrp,
                              i_part,
                              pn_cell[i_part],
                              pn_face[i_part],
                              -1, // pn_edge[i_part],
                              pn_vtx[i_part],
                              pcell_face_idx[i_part],
                              pcell_face[i_part],
                              NULL, //pface_edge_idx[i_part],
                              NULL, //pface_edge[i_part],
                              NULL, //pedge_vtx[i_part],
                              pface_vtx_idx[i_part],
                              pface_vtx[i_part],
                              pcell_ln_to_gn[i_part],
                              pface_ln_to_gn[i_part],
                              NULL, //pedge_ln_to_gn[i_part],
                              pvtx_ln_to_gn[i_part],
                              pvtx_coord[i_part]);

    PDM_extract_part_selected_lnum_set(extrp,
                                       i_part,
                                       pn_select_cell[i_part],
                                       selected_l_num[i_part]);

    // PDM_log_trace_array_int(selected_l_num[i_part], pn_select_cell[i_part], "selected_l_num ::");

  }


  PDM_extract_part_compute(extrp);


  PDM_extract_part_free(extrp);


  for (int i_part = 0; i_part < n_part_zones; i_part++){
    free(cell_center       [i_part]);
    free(selected_g_num    [i_part]);
    free(selected_l_num    [i_part]);
    free(selected_g_num_idx[i_part]);
  }
  free(cell_center);
  free(selected_g_num);
  free(selected_l_num);
  free(selected_g_num_idx);
  free(pn_cell);
  free(pn_face);
  free(pn_vtx);
  free(pn_select_cell);

  free(pcell_ln_to_gn  );
  free(pface_ln_to_gn  );
  free(pvtx_ln_to_gn  );
  free(pcell_face    );
  free(pcell_face_idx);
  free(pface_vtx     );
  free(pface_vtx_idx );
  free(pvtx_coord    );

  PDM_multipart_free(mpart_id);
  PDM_dcube_gen_free(dcube);
  PDM_dmesh_free(dm);

  PDM_MPI_Finalize();

  printf("-- Fin test\n");
  fflush(stdout);

  return 0;
}
