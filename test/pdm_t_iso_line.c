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
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length)
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
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static
inline
double
_unit_circle
(
 double x1,
 double x2
)
{
  return x1 * x1 + x2 * x2 - 0.125;
}

static
inline
void
_unit_circle_gradient
(
 double  x1,
 double  x2,
 double *df_dx1,
 double *df_dx2,
)
{
  *df_dx1 = 2*x1;
  *df_dx2 = 2*x2;
}


static
void
_iso_line
(
 PDM_MPI_Comm  comm,
 int           n_face,
 int           n_edge,
 int           n_vtx,
 int          *pface_edge_idx,
 int          *pface_edge,
 int          *pedge_vtx_idx,
 int          *pedge_vtx,
 PDM_g_num_t  *pedge_ln_to_gn,
 PDM_g_num_t  *pvtx_ln_to_gn,
 double       *pvtx_coord
)
{
  int i_rank;
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);

  PDM_UNUSED(comm);
  PDM_UNUSED(n_face);
  // PDM_UNUSED(n_edge);
  PDM_UNUSED(n_vtx);
  PDM_UNUSED(pface_edge_idx);
  PDM_UNUSED(pface_edge);
  // PDM_UNUSED(pedge_vtx_idx);
  // PDM_UNUSED(pedge_vtx);
  PDM_UNUSED(pedge_ln_to_gn);
  PDM_UNUSED(pvtx_ln_to_gn);
  PDM_UNUSED(pvtx_coord);

  char filename[999];
  sprintf(filename, "out_equi_vtx_coord_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            n_vtx,
                            pvtx_coord,
                            pvtx_ln_to_gn,
                            NULL);

  /*
   *  Tag edges that cross the iso-line,
   *  compute the intersection point
   *  and the gradient at that point
   */
  int    *tag_edge   = PDM_array_zeros_int(n_edge);
  double *edge_coord = (double *) malloc(sizeof(double) * n_edge * 3);
  double *edge_grad  = (double *) malloc(sizeof(double) * n_edge * 3);

  for (int i = 0; i < n_edge; i++) {

    int i_vtx1 = pedge_vtx[2*i  ]-1;
    int i_vtx2 = pedge_vtx[2*i+1]-1;

    double x1 = pvtx_coord[3*i_vtx1  ];
    double y1 = pvtx_coord[3*i_vtx1+1];

    double x2 = pvtx_coord[3*i_vtx2  ];
    double y2 = pvtx_coord[3*i_vtx2+1];

    double val1 = _unit_circle(x1, y1);
    double val2 = _unit_circle(x2, y2);

    int sgn1 = PDM_SIGN(val1);
    int sgn2 = PDM_SIGN(val2);

    if (sgn1 != sgn2) {
      tag_edge[i] = 1;

      double grad1[2], grad2[2];
      _unit_circle_gradient(x1, y1, &grad1[0], &grad1[1]);
      _unit_circle_gradient(x2, y2, &grad2[0], &grad2[1]);

      double t = val1 / (val1 - val2);

      edge_coord[3*i  ] = (1. - t)*x1 + t*y1;
      edge_coord[3*i+1] = (1. - t)*x2 + t*y2;
      edge_coord[3*i+2] = 0.;

      edge_gradient[3*i  ] = (1. - t)*grad1[0] + t*grad1[1];
      edge_gradient[3*i+1] = (1. - t)*grad2[0] + t*grad2[1];
      edge_gradient[3*i+2] = 0.;
    }

  }


  /*
   *
   */

  free(tag_edge);
  free(edge_coord);
  free(edge_gradient);

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

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length);

  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  /*
   *  Create distributed cube
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         -0.5,
                                                         -0.5,
                                                         0.,
                                                         PDM_MESH_NODAL_TRIA3,
                                                         // PDM_MESH_NODAL_QUAD4,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

  if(1 == 1) {
    // PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);

  PDM_dmesh_t* dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);

  int         *dface_edge_idx;
  PDM_g_num_t *dface_edge;
  int dn_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                           &dface_edge,
                                           &dface_edge_idx,
                                           PDM_OWNERSHIP_KEEP);
  int         *dedge_vtx_idx;
  PDM_g_num_t *dedge_vtx;
  int dn_edge  = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                           &dedge_vtx,
                                           &dedge_vtx_idx,
                                           PDM_OWNERSHIP_KEEP);

  if(0 == 1) {
    PDM_log_trace_connectivity_long(dface_edge_idx, dface_edge, dn_face, "dface_edge ::");
    PDM_log_trace_connectivity_long(dedge_vtx_idx , dedge_vtx , dn_edge, "dedge_vtx  ::");
  }

  PDM_UNUSED(dn_face);
  PDM_UNUSED(dn_edge);
  PDM_UNUSED(dn_vtx);
  PDM_UNUSED(dvtx_coord);

  /*
   * Select gnum that contains iso-surface
   */
  PDM_g_num_t* distrib_edge = NULL;
  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_EDGE  , &distrib_edge);
  assert(distrib_edge != NULL);

  PDM_g_num_t* distrib_face = NULL;
  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE  , &distrib_face);
  assert(distrib_face != NULL);

  PDM_g_num_t* edge_ln_to_gn = (PDM_g_num_t * ) malloc( dn_edge * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_edge; ++i) {
    edge_ln_to_gn[i] = distrib_edge[i_rank] + i + 1;
  }

  int          pn_vtx           = 0;
  PDM_g_num_t *pvtx_ln_to_gn  = NULL;
  int         *pedge_vtx_idx    = NULL;
  int         *pedge_vtx        = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_edge,
                                                           dedge_vtx_idx,
                                                           dedge_vtx,
                                                           dn_edge,
                                     (const PDM_g_num_t *) edge_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pedge_vtx_idx,
                                                           &pedge_vtx);

  double **tmp_pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        vtx_distrib,
                                        dvtx_coord,
                                        &pn_vtx,
                 (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coord);
  double* pvtx_coord = tmp_pvtx_coord[0];
  free(tmp_pvtx_coord);

  /*
   * Select edge
   */
  int *dedge_tag = (int * ) malloc(dn_edge * sizeof(int));

  for(int i = 0; i < dn_edge; ++i) {

    int i_vtx1 = pedge_vtx[2*i  ]-1;
    int i_vtx2 = pedge_vtx[2*i+1]-1;

    dedge_tag[i] = 0;

    double x1 = pvtx_coord[3*i_vtx1  ];
    double y1 = pvtx_coord[3*i_vtx1+1];
    // int x3 = pvtx_coord[3*i_vtx1+2];

    double x2 = pvtx_coord[3*i_vtx2  ];
    double y2 = pvtx_coord[3*i_vtx2+1];

    double val1 = _unit_circle(x1, y1);
    double val2 = _unit_circle(x2, y2);

    int sgn1 = PDM_SIGN(val1);
    int sgn2 = PDM_SIGN(val2);

    if(sgn1 * sgn2 < 0) {
      dedge_tag[i] = 1;
    }
  }

  int* val = malloc(pn_vtx * sizeof(int));
  printf("pn_vtx = %i \n", pn_vtx);
  for(int i = 0; i < pn_vtx; ++i) {
    double x1 = pvtx_coord[3*i  ];
    double y1 = pvtx_coord[3*i+1];
    val[i]    = PDM_SIGN( _unit_circle(x1, y1) );
  }

  char filename[999];
  sprintf(filename, "out_edge_tag_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            pn_vtx,
                            pvtx_coord,
                            NULL,
                            val);
  free(val);


  // PDM_log_trace_array_int(dedge_tag, dn_edge, "dedge_tag");
  free(pvtx_coord);

  /*
   * block_to_part on dface_edge
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_edge,
                               (const PDM_g_num_t **) &dface_edge,
                                                      &dface_edge_idx[dn_face],
                                                      1,
                                                      comm);

  int strid_one = 1;
  int **tmp_dface_edge_tag = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST_INTERLACED,
                         &strid_one,
            (void *  )   dedge_tag,
            (int  ***)   NULL,
            (void ***)  &tmp_dface_edge_tag);
  int *dface_edge_tag = tmp_dface_edge_tag[0];
  free(tmp_dface_edge_tag);
  free(dedge_tag);

  int         *dface_tag        = malloc(dn_face * sizeof(int        ));
  PDM_g_num_t *face_to_extract_gnum = malloc(dn_face * sizeof(PDM_g_num_t));
  int  n_face_tag = 0;
  for(int i = 0; i < dn_face; ++i) {
    dface_tag[i] = 0;
    for(int idx_face = dface_edge_idx[i]; idx_face < dface_edge_idx[i+1]; ++idx_face) {
      if(dface_edge_tag[idx_face] == 1) {
        dface_tag[i] = 1;
        face_to_extract_gnum[n_face_tag++] = distrib_face[i_rank] + i + 1;
        break;
      }
    }
  }

  if(0 == 1) {
    PDM_log_trace_array_int (dedge_tag, dn_face, "dface_tag");
    PDM_log_trace_array_long(face_to_extract_gnum, n_face_tag, "face_to_extract_gnum");
  }

  PDM_block_to_part_free(btp);

  free(edge_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(pedge_vtx_idx);
  free(pedge_vtx);
  free(dface_edge_tag);

  /*
   *  Visu
   */
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;

  for(int i = 0; i < dface_edge_idx[dn_face]; ++i) {
    dface_edge[i] = PDM_ABS(dface_edge[i]);
  }

  PDM_deduce_combine_connectivity(comm,
                                  distrib_face,
                                  distrib_edge,
                                  dface_edge_idx,
                                  dface_edge,
                                  dedge_vtx_idx,
                                  dedge_vtx,
                                  1,
                                  &dface_vtx_idx,
                                  &dface_vtx);

  // PDM_log_trace_connectivity_long(dface_vtx_idx, dface_vtx, dn_face, "dface_edge ::");
  for(int i = 0; i < dface_vtx_idx[dn_face]; ++i) {
    dface_vtx[i] = PDM_ABS(dface_vtx[i]);
  }

  // PDM_log_trace_array_int(dface_vtx_idx, dn_face+1,"dface_vtx_idx :: " );
  // PDM_log_trace_array_int(dface_edge_idx, dn_face+1,"dface_edge_idx :: " );

  int          pn_extract_vtx           = 0;
  PDM_g_num_t *pn_extract_vtx_ln_to_gn  = NULL;
  int         *pface_vtx_idx    = NULL;
  int         *pface_vtx        = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_vtx_idx,
                                                           dface_vtx,
                                                           n_face_tag,
                                     (const PDM_g_num_t *) face_to_extract_gnum,
                                                           &pn_extract_vtx,
                                                           &pn_extract_vtx_ln_to_gn,
                                                           &pface_vtx_idx,
                                                           &pface_vtx);

  double **tmp_pvtx_extract_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        vtx_distrib,
                                        dvtx_coord,
                                        &pn_extract_vtx,
                 (const PDM_g_num_t **) &pn_extract_vtx_ln_to_gn,
                                        &tmp_pvtx_extract_coord);
  double* pvtx_extract_coord = tmp_pvtx_extract_coord[0];
  free(tmp_pvtx_extract_coord);


  sprintf(filename, "out_mesh_%2.2d.vtk", i_rank);
  PDM_vtk_write_polydata(filename,
                         pn_extract_vtx,
                         pvtx_extract_coord,
                         pn_extract_vtx_ln_to_gn,//NULL,
                         n_face_tag,
                         pface_vtx_idx,
                         pface_vtx,
                         face_to_extract_gnum,
                         NULL);



  free(pn_extract_vtx_ln_to_gn);



  double* face_center = (double *) malloc( 3 * n_face_tag * sizeof(double));
  for(int i_face = 0; i_face < n_face_tag; ++i_face) {
    face_center[3*i_face  ] = 0.;
    face_center[3*i_face+1] = 0.;
    face_center[3*i_face+2] = 0.;
    int nvtx_on_face = pface_vtx_idx[i_face+1] - pface_vtx_idx[i_face];
    for(int idx_vtx = pface_vtx_idx[i_face]; idx_vtx < pface_vtx_idx[i_face+1]; ++idx_vtx) {
      int i_vtx = pface_vtx[idx_vtx]-1;
      face_center[3*i_face  ] += pvtx_extract_coord[3*i_vtx  ];
      face_center[3*i_face+1] += pvtx_extract_coord[3*i_vtx+1];
      face_center[3*i_face+2] += pvtx_extract_coord[3*i_vtx+2];
    }
    double inv = 1./(double) nvtx_on_face;
    face_center[3*i_face  ] = face_center[3*i_face  ] * inv;
    face_center[3*i_face+1] = face_center[3*i_face+1] * inv;
    face_center[3*i_face+2] = face_center[3*i_face+2] * inv;
  }

  sprintf(filename, "out_equi_face_coord_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            n_face_tag,
                            face_center,
                            NULL,
                            NULL);


  free(pface_vtx_idx);
  free(pface_vtx);


  free(pvtx_extract_coord);

  /*
   * Rebuild partition that contains faces and reequilibrate
   */
  PDM_gen_gnum_t* gnum_equi = PDM_gnum_create(3, 1, PDM_FALSE, 0., comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_coords(gnum_equi, 0, n_face_tag, face_center, NULL);
  PDM_gnum_compute(gnum_equi);
  PDM_g_num_t* child_equi_face_gnum = PDM_gnum_get(gnum_equi, 0);
  PDM_gnum_free(gnum_equi);
  free(face_center);

  /*
   * Equilibrage avec le part_to_block
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                       1.,
                                                       &child_equi_face_gnum,
                                                       NULL,
                                                       &n_face_tag,
                                                       1,
                                                       comm);

  int n_face_equi = PDM_part_to_block_n_elt_block_get (ptb);
  PDM_g_num_t *block_face_g_num_child_equi = PDM_part_to_block_block_gnum_get (ptb);

  PDM_g_num_t *block_face_equi_parent_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
               (void **) &face_to_extract_gnum,
                          NULL,
               (void **) &block_face_equi_parent_g_num);

  if(0 == 1) {
    PDM_log_trace_array_long(block_face_equi_parent_g_num, n_face_equi, "block_face_equi_parent_g_num ::");
    PDM_log_trace_array_long(block_face_g_num_child_equi , n_face_equi, "block_face_g_num_child_equi  ::");
  }

  /*
   * Je prepare tout pour mon petit Bastien
   */
  int          pn_edge_equi        = 0;
  PDM_g_num_t *pequi_parent_edge_ln_to_gn = NULL;
  int         *pequi_face_edge_idx = NULL;
  int         *pequi_face_edge     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_edge_idx,
                                                           dface_edge,
                                                           n_face_equi,
                                     (const PDM_g_num_t *) block_face_equi_parent_g_num,
                                                           &pn_edge_equi,
                                                           &pequi_parent_edge_ln_to_gn,
                                                           &pequi_face_edge_idx,
                                                           &pequi_face_edge);

  PDM_gen_gnum_t* gnum_edge = PDM_gnum_create(3, 1, PDM_FALSE, 0., comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_parents(gnum_edge, 0, pn_edge_equi, pequi_parent_edge_ln_to_gn);
  PDM_gnum_compute(gnum_edge);
  PDM_g_num_t* pequi_edge_ln_to_gn = PDM_gnum_get(gnum_edge, 0);
  PDM_gnum_free(gnum_edge);


  int          pn_vtx_equi        = 0;
  PDM_g_num_t *pequi_parent_vtx_ln_to_gn = NULL;
  int         *pequi_edge_vtx_idx = NULL;
  int         *pequi_edge_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_edge,
                                                           dedge_vtx_idx,
                                                           dedge_vtx,
                                                           pn_edge_equi,
                                     (const PDM_g_num_t *) pequi_parent_edge_ln_to_gn,
                                                           &pn_vtx_equi,
                                                           &pequi_parent_vtx_ln_to_gn,
                                                           &pequi_edge_vtx_idx,
                                                           &pequi_edge_vtx);

  PDM_gen_gnum_t* gnum_vtx = PDM_gnum_create(3, 1, PDM_FALSE, 0., comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_parents(gnum_vtx, 0, pn_vtx_equi, pequi_parent_vtx_ln_to_gn);
  PDM_gnum_compute(gnum_vtx);
  PDM_g_num_t* pequi_vtx_ln_to_gn = PDM_gnum_get(gnum_vtx, 0);
  PDM_gnum_free(gnum_vtx);

  double **tmp_pequi_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        vtx_distrib,
                                        dvtx_coord,
                                        &pn_vtx_equi,
                 (const PDM_g_num_t **) &pequi_parent_vtx_ln_to_gn,
                                        &tmp_pequi_vtx_coord);
  double* pequi_vtx_coord = tmp_pequi_vtx_coord[0];
  free(tmp_pequi_vtx_coord);

  _iso_line(comm,
            n_face_equi,
            pn_edge_equi,
            pn_vtx_equi,
            pequi_face_edge_idx,
            pequi_face_edge,
            pequi_edge_vtx_idx,
            pequi_edge_vtx,
            pequi_edge_ln_to_gn,
            pequi_vtx_ln_to_gn,
            pequi_vtx_coord);

  free(pequi_edge_ln_to_gn);
  free(pequi_vtx_ln_to_gn);


  free(pequi_vtx_coord);

  free(pequi_parent_vtx_ln_to_gn);
  free(pequi_edge_vtx_idx);
  free(pequi_edge_vtx);

  free(pequi_parent_edge_ln_to_gn);
  free(pequi_face_edge_idx);
  free(pequi_face_edge);


  free(block_face_equi_parent_g_num);


  PDM_part_to_block_free(ptb);
  free(dface_vtx_idx);
  free(dface_vtx);
  free(child_equi_face_gnum);


  free(face_to_extract_gnum);
  free(dface_tag);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_dcube_nodal_gen_free(dcube);


  PDM_MPI_Finalize();

  return 0;
}
