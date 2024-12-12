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
#include "pdm_part_mesh.h"

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
           PDM_g_num_t   *n_vtx_seg,
           double        *length,
           int           *n_part,
           int           *post,
           int           *part_method,
           int           *local)
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
    else if (strcmp(argv[i], "-local") == 0) {
      *local = 1;
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

  PDM_g_num_t       n_vtx_seg   = 10;
  double            length      = 1.;
  int               n_part      = 1;
  int               post        = 0;
  PDM_split_dual_t  part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  int               local       = 0;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
     (int *) &part_method,
             &local);

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
  PDM_dmesh_t* dm = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                     dn_cell,
                                     dn_face,
                                     0, // dn_edge
                                     dn_vtx,
                                     comm);

  PDM_dmesh_vtx_coord_set(dm,
                          dvtx_coord,
                          PDM_OWNERSHIP_USER);


  PDM_dmesh_connectivity_set(dm,
                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
                             dface_vtx,
                             dface_vtx_idx,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dm,
                             PDM_CONNECTIVITY_TYPE_FACE_CELL,
                             dface_cell,
                             NULL,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_bound_set(dm,
                      PDM_BOUND_TYPE_FACE,
                      n_face_group,
                      dface_group,
                      dface_group_idx,
                      PDM_OWNERSHIP_USER);

  /*
   * Partitionnement
   */
  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_domains,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart, -1, "PDM_PART_RENUM_CELL_NONE", NULL, "PDM_PART_RENUM_FACE_NONE");
  PDM_multipart_dmesh_set(mpart, 0, dm);
  PDM_multipart_compute(mpart);

  /*
   * Get the partition domain
   */
  int i_domain = 0;

  double      **cell_center;
  int         **selected_l_num;
  PDM_g_num_t **pcell_ln_to_gn;
  PDM_g_num_t **pface_ln_to_gn;
  PDM_g_num_t **pvtx_ln_to_gn;
  int          *pn_cell;
  int          *pn_face;
  int          *pn_vtx;
  int          *pn_select_cell;
  int         **pcell_face;
  int         **pcell_face_idx;
  int         **pface_vtx;
  int         **pface_vtx_idx;
  double      **pvtx_coord;
  PDM_malloc(selected_l_num, n_part_domains, int         *);
  PDM_malloc(pcell_ln_to_gn, n_part_domains, PDM_g_num_t *);
  PDM_malloc(pface_ln_to_gn, n_part_domains, PDM_g_num_t *);
  PDM_malloc(pvtx_ln_to_gn,  n_part_domains, PDM_g_num_t *);
  PDM_malloc(pn_cell,        n_part_domains, int          );
  PDM_malloc(pn_face,        n_part_domains, int          );
  PDM_malloc(pn_vtx,         n_part_domains, int          );
  PDM_malloc(pn_select_cell, n_part_domains, int          );
  PDM_malloc(pcell_face,     n_part_domains, int         *);
  PDM_malloc(pcell_face_idx, n_part_domains, int         *);
  PDM_malloc(pface_vtx,      n_part_domains, int         *);
  PDM_malloc(pface_vtx_idx,  n_part_domains, int         *);
  PDM_malloc(pvtx_coord,     n_part_domains, double      *);

  for (int i_part = 0; i_part < n_part_domains; i_part++){

    PDM_g_num_t* cell_ln_to_gn = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    i_domain,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    &cell_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    int *cell_face     = NULL;
    int *cell_face_idx = NULL;
    int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                     i_domain,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                     &cell_face_idx,
                                                     &cell_face,
                                                     PDM_OWNERSHIP_KEEP);

    int *face_vtx     = NULL;
    int *face_vtx_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        i_domain,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                        &face_vtx_idx,
                                        &face_vtx,
                                        PDM_OWNERSHIP_KEEP);

    PDM_g_num_t* face_ln_to_gn = NULL;
    int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                                 i_domain,
                                                 i_part,
                                                 PDM_MESH_ENTITY_FACE,
                                                 &face_ln_to_gn,
                                                 PDM_OWNERSHIP_KEEP);

    PDM_g_num_t* vtx_ln_to_gn = NULL;
    int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                i_domain,
                                                i_part,
                                                PDM_MESH_ENTITY_VTX,
                                                &vtx_ln_to_gn,
                                                PDM_OWNERSHIP_KEEP);

    double *vtx = NULL;
    PDM_multipart_part_vtx_coord_get(mpart,
                                     i_domain,
                                     i_part,
                                     &vtx,
                                     PDM_OWNERSHIP_KEEP);


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

    // char filename[999];
    // sprintf(filename, "mesh_before_extract_%3.3d_%3.3d.vtk", i_part, i_rank);
    // PDM_vtk_write_polydata(filename,
    //                        pn_vtx[i_part],
    //                        pvtx_coord[i_part],
    //                        pvtx_ln_to_gn[i_part],
    //                        pn_face[i_part],
    //                        pface_vtx_idx[i_part],
    //                        pface_vtx[i_part],
    //                        pface_ln_to_gn[i_part],
    //                        NULL);
  }

  /*
   * Compute cell centers
   */
  PDM_part_geom_cell_center(n_part,
                            pn_cell,
                            NULL,
                            pcell_face_idx,
                            pcell_face,
                            NULL,
                            NULL,
                            pface_vtx_idx,
                            pface_vtx,
                            NULL,
                            pvtx_coord,
                            &cell_center);
    // double *face_center;
    // PDM_malloc(face_center, 3 * n_face ,double);

    // for(int i_face = 0; i_face < n_face; ++i_face) {
    //   face_center[3*i_face  ] = 0.;
    //   face_center[3*i_face+1] = 0.;
    //   face_center[3*i_face+2] = 0.;
    //   int n_vtx_on_face = face_vtx_idx[i_face+1] - face_vtx_idx[i_face];
    //   for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
    //     int i_vtx = PDM_ABS(face_vtx[idx_vtx])-1;
    //     face_center[3*i_face  ] += vtx[3*i_vtx  ];
    //     face_center[3*i_face+1] += vtx[3*i_vtx+1];
    //     face_center[3*i_face+2] += vtx[3*i_vtx+2];
    //   }
    //   face_center[3*i_face  ] = face_center[3*i_face  ] / n_vtx_on_face;
    //   face_center[3*i_face+1] = face_center[3*i_face+1] / n_vtx_on_face;
    //   face_center[3*i_face+2] = face_center[3*i_face+2] / n_vtx_on_face;
    // }

    // PDM_malloc(cell_center[i_part], 3 * n_cell ,double);

    // for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

    //   cell_center[i_part][3*i_cell  ] = 0.;
    //   cell_center[i_part][3*i_cell+1] = 0.;
    //   cell_center[i_part][3*i_cell+2] = 0.;

    //   int n_face_on_cell = cell_face_idx[i_cell+1] - cell_face_idx[i_cell];

    //   for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {

    //     int i_face = PDM_ABS(cell_face[idx_face])-1;
    //     cell_center[i_part][3*i_cell  ] += face_center[3*i_face  ];
    //     cell_center[i_part][3*i_cell+1] += face_center[3*i_face+1];
    //     cell_center[i_part][3*i_cell+2] += face_center[3*i_face+2];
    //   }
    //   cell_center[i_part][3*i_cell  ] = cell_center[i_part][3*i_cell  ] / n_face_on_cell;
    //   cell_center[i_part][3*i_cell+1] = cell_center[i_part][3*i_cell+1] / n_face_on_cell;
    //   cell_center[i_part][3*i_cell+2] = cell_center[i_part][3*i_cell+2] / n_face_on_cell;
    // }

    // PDM_free(face_center);

  /**
   *  Extract cells intersecting a given box
   */
  double bbox[6];
  // bbox[0] = 0.3;
  // bbox[1] = 0.3;
  // bbox[2] = 0.35;
  // bbox[3] = 0.7;
  // bbox[4] = 0.7;
  // bbox[5] = 0.65;
  bbox[0] = 0.5;
  bbox[1] = 0.65;
  bbox[2] = 0.75;
  bbox[3] = 1.25;
  bbox[4] = 1.25;
  bbox[5] = 1.25;

  for (int i_part = 0; i_part < n_part_domains; i_part++){
    PDM_malloc(selected_l_num[i_part], pn_cell[i_part], int);

    /*
     * Sub-part
     */
    int n_select_cell = 0;
    for(int i_cell = 0; i_cell < pn_cell[i_part]; ++i_cell) {

      int inside = 1;
      for(int i = 0; i < 3; ++i) {
        if (cell_center[i_part][3*i_cell+i] > bbox[i+3] || cell_center[i_part][3*i_cell+i] < bbox[i]) {
          inside = 0;
        }
      }
      if(inside == 1) {
        selected_l_num[i_part][n_select_cell] = i_cell+1;
        n_select_cell++;
      }

      // if(cell_ln_to_gn[i_cell] == 5) {
      // if(cell_ln_to_gn[i_cell]%2 == 1) {
      //   selected_l_num[i_part][n_select_cell]     = i_cell;
      //   n_select_cell++;
      // }


    }

    PDM_realloc(selected_l_num[i_part], selected_l_num[i_part], n_select_cell, int);
    pn_select_cell[i_part] = n_select_cell;

  }
  //PDM_free(dface_join_idx);

  /*
   * Extract
   */
  int n_part_out       = 1;
  int n_bound          = 0;
  int n_fake_group_vtx = 2;

  PDM_extract_part_kind_t extract_kind = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
  if (local) {
    extract_kind = PDM_EXTRACT_PART_KIND_LOCAL;
  }
  // PDM_split_dual_t        split_dual_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  PDM_split_dual_t        split_dual_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  // PDM_split_dual_t        split_dual_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  PDM_extract_part_t* extrp = PDM_extract_part_create(3,
                                                      n_part,
                                                      n_part_out,
                                                      extract_kind,
                                                      split_dual_method,
                                                      PDM_TRUE, // compute_child_gnum
                                                      PDM_OWNERSHIP_KEEP,
                                                      comm);

  int         **fake_group_vtx;
  int         **fake_group_face;
  int         **fake_group_face_idx;
  PDM_g_num_t **fake_face_group_ln_to_gn;
  PDM_malloc(fake_group_vtx,           n_part, int         *);
  PDM_malloc(fake_group_face,          n_part, int         *);
  PDM_malloc(fake_group_face_idx,      n_part, int         *);
  PDM_malloc(fake_face_group_ln_to_gn, n_part, PDM_g_num_t *);

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
                                       selected_l_num[i_part],
                                       PDM_OWNERSHIP_USER);


    int         *group_face_idx      = NULL;
    int         *group_face          = NULL;
    PDM_g_num_t *face_group_ln_to_gn = NULL;
    PDM_multipart_group_get(mpart,
                            i_domain,
                            i_part,
                            PDM_MESH_ENTITY_FACE,
                            &n_bound,
                            &group_face_idx,
                            &group_face,
                            &face_group_ln_to_gn,
                            PDM_OWNERSHIP_KEEP);

    /** Vertex groups **/

    /* Add 2 fake groups */
    PDM_malloc(fake_group_vtx[i_part], pn_vtx[i_part], int);

    for (int i_vtx = 0; i_vtx < pn_vtx[i_part]; ++i_vtx){
      fake_group_vtx[i_part][i_vtx] = i_vtx+1;
    }

    PDM_extract_part_n_group_set(extrp,
                                 PDM_BOUND_TYPE_VTX,
                                 n_fake_group_vtx);

    for(int i_group = 0; i_group < n_fake_group_vtx; ++i_group) {
      PDM_extract_part_part_group_set(extrp,
                                      i_part,
                                      i_group,
                                      PDM_BOUND_TYPE_VTX,
                                      pn_vtx[i_part],
                                      fake_group_vtx[i_part],
                                      pvtx_ln_to_gn[i_part]);
    }

    /** Face groups **/

    /* Add n_bound+1 groups */
    int fake_group_n_face = group_face_idx[n_bound]-group_face_idx[n_bound-1];

    PDM_malloc(fake_group_face[i_part],          group_face_idx[n_bound]+fake_group_n_face, int        );
    PDM_malloc(fake_group_face_idx[i_part],      n_bound+2                                , int        );
    PDM_malloc(fake_face_group_ln_to_gn[i_part], group_face_idx[n_bound]+fake_group_n_face, PDM_g_num_t);

    fake_group_face_idx[i_part][0] = 0;

    /* by copying existing groups */
    for(int i_group = 0; i_group < n_bound; ++i_group) {
      fake_group_face_idx[i_part][i_group+1] = group_face_idx[i_group+1];
      for (int i_face = group_face_idx[i_group]; i_face < group_face_idx[i_group+1]; ++i_face){
        fake_group_face[i_part][i_face]          = group_face[i_face];
        fake_face_group_ln_to_gn[i_part][i_face] = face_group_ln_to_gn[i_face];
      }
    }

    /* and duplicating the last group */
    fake_group_face_idx[i_part][n_bound+1] = fake_group_face_idx[i_part][n_bound] + fake_group_n_face;
    for (int i_face = fake_group_face_idx[i_part][n_bound]; i_face < fake_group_face_idx[i_part][n_bound+1]; ++i_face){
      int j_face = i_face - fake_group_n_face;
      fake_group_face[i_part][i_face]          = group_face[j_face];
      fake_face_group_ln_to_gn[i_part][i_face] = face_group_ln_to_gn[j_face];
    }

    PDM_extract_part_n_group_set(extrp,
                                 PDM_BOUND_TYPE_FACE,
                                 n_bound+1);

    for(int i_group = 0; i_group < n_bound+1; ++i_group) {
      PDM_extract_part_part_group_set(extrp,
                                      i_part,
                                      i_group,
                                      PDM_BOUND_TYPE_FACE,
                                      fake_group_face_idx[i_part][i_group+1]-fake_group_face_idx[i_part][i_group],
                                      &fake_group_face[i_part][group_face_idx[i_group]],
                                      &fake_face_group_ln_to_gn[i_part][group_face_idx[i_group]]);
    }

    // PDM_log_trace_array_int(selected_l_num[i_part], pn_select_cell[i_part], "selected_l_num ::");

  }


  PDM_extract_part_compute(extrp);

  int          *pn_extract_face;
  int          *pn_extract_vtx;
  int         **pextract_face_vtx;
  int         **pextract_face_vtx_idx;
  double      **pextract_vtx;
  PDM_g_num_t **pextract_face_ln_to_gn;
  int         **pextract_face_group;
  PDM_g_num_t **pextract_vtx_ln_to_gn;
  PDM_malloc(pn_extract_face,        n_part_out, int          );
  PDM_malloc(pn_extract_vtx,         n_part_out, int          );
  PDM_malloc(pextract_face_vtx,      n_part_out, int         *);
  PDM_malloc(pextract_face_vtx_idx,  n_part_out, int         *);
  PDM_malloc(pextract_vtx,           n_part_out, double      *);
  PDM_malloc(pextract_face_ln_to_gn, n_part_out, PDM_g_num_t *);
  PDM_malloc(pextract_face_group,    n_part_out, int         *);
  PDM_malloc(pextract_vtx_ln_to_gn,  n_part_out, PDM_g_num_t *);


  for(int i_part = 0; i_part < n_part_out; ++i_part) {

    pn_extract_face[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                            i_part,
                                                            PDM_MESH_ENTITY_FACE);

    pn_extract_vtx[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                           i_part,
                                                           PDM_MESH_ENTITY_VTX);

    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                      &pextract_face_vtx[i_part],
                                      &pextract_face_vtx_idx[i_part],
                                      PDM_OWNERSHIP_KEEP);

    PDM_extract_part_vtx_coord_get(extrp,
                                   i_part,
                                   &pextract_vtx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  &pextract_face_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_VTX,
                                  &pextract_vtx_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    PDM_malloc(pextract_face_group[i_part], pn_extract_face[i_part], int);
    for(int i_face = 0; i_face < pn_extract_face[i_part]; ++i_face) {
      pextract_face_group[i_part][i_face] = -1;
    }

    for(int i_group = 0; i_group < n_fake_group_vtx; ++i_group) {
      int          pn_extract_group_entity               = 0;
      int         *pextract_group_entity                 = NULL;
      PDM_g_num_t *pextract_group_entity_ln_to_gn        = NULL;
      PDM_g_num_t *pextract_group_entity_parent_ln_to_gn = NULL;
      PDM_extract_part_group_get(extrp,
                                 PDM_BOUND_TYPE_VTX,
                                 i_part,
                                 i_group,
                                 &pn_extract_group_entity,
                                 &pextract_group_entity,
                                 &pextract_group_entity_ln_to_gn,
                                 &pextract_group_entity_parent_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

      /* Check extraction */
      if (i_group == 0){
        assert(pn_extract_group_entity == 72);
      }
      if (i_group == 1){
        assert(pn_extract_group_entity == 72);
      }

    }

    for(int i_group = 0; i_group < n_bound+1; ++i_group) {
      int          pn_extract_group_entity               = 0;
      int         *pextract_group_entity                 = NULL;
      PDM_g_num_t *pextract_group_entity_ln_to_gn        = NULL;
      PDM_g_num_t *pextract_group_entity_parent_ln_to_gn = NULL;
      PDM_extract_part_group_get(extrp,
                                 PDM_BOUND_TYPE_FACE,
                                 i_part,
                                 i_group,
                                 &pn_extract_group_entity,
                                 &pextract_group_entity,
                                 &pextract_group_entity_ln_to_gn,
                                 &pextract_group_entity_parent_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

      for(int idx_entity = 0; idx_entity < pn_extract_group_entity; ++idx_entity) {
        int i_face = pextract_group_entity[idx_entity]-1;
        pextract_face_group[i_part][i_face] = i_group;
      }

      /* Check extraction */
      if (i_group == n_bound-1){
        assert(pn_extract_group_entity == 10);
      }
      if (i_group == n_bound){
        assert(pn_extract_group_entity == 10);
      }

    }

    // PDM_g_num_t *pextract_parent_cell_ln_to_gn = NULL;
    // PDM_g_num_t *pextract_cell_ln_to_gn = NULL;
    // int n_cell = PDM_extract_part_ln_to_gn_get(extrp,
    //                                            i_part,
    //                                            PDM_MESH_ENTITY_CELL,
    //                                            &pextract_cell_ln_to_gn,
    //                                            PDM_OWNERSHIP_KEEP);

    // n_cell = PDM_extract_part_parent_ln_to_gn_get(extrp,
    //                                                   i_part,
    //                                                   PDM_MESH_ENTITY_CELL,
    //                                                   &pextract_parent_cell_ln_to_gn,
    //                                                   PDM_OWNERSHIP_KEEP);

    // PDM_log_trace_array_long(pextract_cell_ln_to_gn, n_cell, "pextract_cell_ln_to_gn ::");
    // PDM_log_trace_array_long(pextract_parent_cell_ln_to_gn, n_cell, "pextract_parent_cell_ln_to_gn ::");


  }

  /*
   * Export vtk en légende
   */
  if(post) {
    for(int i_part = 0; i_part < n_part_out; ++i_part) {

      char filename[999];
      sprintf(filename, "extract_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_point_cloud(filename,
                                pn_extract_vtx[i_part],
                                pextract_vtx[i_part],
                                NULL, NULL);

      PDM_log_trace_connectivity_int(pextract_face_vtx_idx[i_part],
                                     pextract_face_vtx    [i_part],
                                     pn_extract_face[i_part], " pextract_face_vtx :: ");

      sprintf(filename, "extract_face_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_polydata(filename,
                             pn_extract_vtx        [i_part],
                             pextract_vtx          [i_part],
                             pextract_vtx_ln_to_gn [i_part],
                             pn_extract_face       [i_part],
                             pextract_face_vtx_idx [i_part],
                             pextract_face_vtx     [i_part],
                             pextract_face_ln_to_gn[i_part],
                             pextract_face_group   [i_part]);
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(fake_group_vtx[i_part]);
    PDM_free(fake_group_face[i_part]);
    PDM_free(fake_group_face_idx[i_part]);
    PDM_free(fake_face_group_ln_to_gn[i_part]);
  }
  PDM_free(fake_group_vtx);
  PDM_free(fake_group_face);
  PDM_free(fake_group_face_idx);
  PDM_free(fake_face_group_ln_to_gn);

  PDM_free(pn_extract_face);
  PDM_free(pn_extract_vtx);
  PDM_free(pextract_face_vtx     );
  PDM_free(pextract_face_vtx_idx );
  PDM_free(pextract_vtx          );
  PDM_free(pextract_face_ln_to_gn);
  PDM_free(pextract_vtx_ln_to_gn );

  // Test output as part_mesh
  PDM_part_mesh_t *pmesh = NULL;
  PDM_bool_t pmesh_takes_ownership = PDM_FALSE;
  PDM_extract_part_part_mesh_get(extrp,
                                 &pmesh,
                                 pmesh_takes_ownership);

  PDM_part_mesh_free(pmesh);

  pmesh_takes_ownership = PDM_TRUE;
  PDM_extract_part_part_mesh_get(extrp,
                                         &pmesh,
                                         pmesh_takes_ownership);

  PDM_part_mesh_free(pmesh);

  PDM_extract_part_free(extrp);


  for (int i_part = 0; i_part < n_part_domains; i_part++){
    PDM_free(cell_center        [i_part]);
    PDM_free(selected_l_num     [i_part]);
    PDM_free(pextract_face_group[i_part]);
  }
  PDM_free(pextract_face_group);
  PDM_free(cell_center);
  PDM_free(selected_l_num);
  PDM_free(pn_cell);
  PDM_free(pn_face);
  PDM_free(pn_vtx);
  PDM_free(pn_select_cell);

  PDM_free(pcell_ln_to_gn);
  PDM_free(pface_ln_to_gn);
  PDM_free(pvtx_ln_to_gn );
  PDM_free(pcell_face    );
  PDM_free(pcell_face_idx);
  PDM_free(pface_vtx     );
  PDM_free(pface_vtx_idx );
  PDM_free(pvtx_coord    );

  PDM_multipart_free(mpart);
  PDM_dcube_gen_free(dcube);
  PDM_dmesh_free(dm);


  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize();

  return 0;
}
