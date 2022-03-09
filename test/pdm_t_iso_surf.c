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
#include "pdm_array.h"

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
_unit_circle_gradient
(
 double  x1,
 double  x2,
 double  x3,
 double *df_dx1,
 double *df_dx2,
 double *df_dx3
)
{
  *df_dx1 = 2*x1;
  *df_dx2 = 2*x2;
  *df_dx3 = 2*x3;
}


static void
_dump_vectors
(
 const char   *filename,
 const int     n_pts,
 const double  pts_coord[],
 const double  vector[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", pts_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", n_pts, 2*n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1\n");
  }

  fprintf(f, "POINT_DATA %d\n", n_pts);
  fprintf(f, "VECTORS vector double\n");
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vector[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fclose(f);
}


static
void
_iso_surf
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
  PDM_MPI_Comm_rank(comm, &i_rank);

  // PDM_UNUSED(comm);
  PDM_UNUSED(n_face);
  // PDM_UNUSED(n_edge);
  PDM_UNUSED(n_vtx);
  PDM_UNUSED(pface_edge_idx);
  PDM_UNUSED(pface_edge);
  PDM_UNUSED(pedge_vtx_idx);
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
  return;
  /*
   *  Tag edges that cross the iso-line,
   *  compute the intersection point
   *  and the normal at that point
   */
  int    *tag_edge    = PDM_array_zeros_int(n_edge);
  double *edge_coord  = (double *) malloc(sizeof(double) * n_edge * 3);
  double *edge_normal = (double *) malloc(sizeof(double) * n_edge * 3);

  for (int i = 0; i < n_edge; i++) {

    int i_vtx1 = pedge_vtx[2*i  ]-1;
    int i_vtx2 = pedge_vtx[2*i+1]-1;

    double x1 = pvtx_coord[3*i_vtx1  ];
    double y1 = pvtx_coord[3*i_vtx1+1];
    double z1 = pvtx_coord[3*i_vtx1+2];

    double x2 = pvtx_coord[3*i_vtx2  ];
    double y2 = pvtx_coord[3*i_vtx2+1];
    double z2 = pvtx_coord[3*i_vtx2+2];

    double val1 = _unit_sphere(x1, y1, z1);
    double val2 = _unit_sphere(x2, y2, z2);

    int sgn1 = PDM_SIGN(val1);
    int sgn2 = PDM_SIGN(val2);

    if (sgn1 != sgn2) {
      tag_edge[i] = 1;

      double grad1[3], grad2[3], grad3[3];
      _unit_circle_gradient(x1, y1, z1, &grad1[0], &grad1[1], &grad1[2]);
      _unit_circle_gradient(x2, y2, z2, &grad2[0], &grad2[1], &grad2[2]);

      double t = val1 / (val1 - val2);

      edge_coord[3*i  ] = (1. - t)*x1 + t*x2;
      edge_coord[3*i+1] = (1. - t)*y1 + t*y2;
      edge_coord[3*i+2] = (1. - t)*z1 + t*z2;

      abort();
      double gx = (1. - t)*grad1[0] + t*grad2[0];
      double gy = (1. - t)*grad1[1] + t*grad2[1];
      double gz = (1. - t)*grad1[1] + t*grad2[1];
      double imag = 1. / sqrt(gx*gx + gy*gy);

      edge_normal[3*i  ] = gx * imag;
      edge_normal[3*i+1] = gy * imag;
      edge_normal[3*i+2] = 0.;
    }

    else {
      edge_coord[3*i  ] = 0.;
      edge_coord[3*i+1] = 0.;
      edge_coord[3*i+2] = 0.;

      edge_normal[3*i  ] = 0.;
      edge_normal[3*i+1] = 0.;
      edge_normal[3*i+2] = 0.;
    }

  }


  sprintf(filename, "edge_intersection_%2.2d.vtk", i_rank);
  _dump_vectors (filename,
                 n_edge,
                 edge_coord,
                 edge_normal);


  int n_face_edge_max = 0;
  for (int i = 0; i < n_face; i++) {
    n_face_edge_max = PDM_MAX(n_face_edge_max,
                              pface_edge_idx[i+1] - pface_edge_idx[i]);
  }

  double *isoline_vtx_coord = (double *) malloc(sizeof(double) * n_face * 3);
  for (int i = 0; i < n_face; i++) {

    int n_face_edge = pface_edge_idx[i+1] - pface_edge_idx[i];

  }


  free(tag_edge);
  free(edge_coord);
  free(edge_normal);
  free(isoline_vtx_coord);

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
                                                         -0.5,
                                                         PDM_MESH_NODAL_HEXA8,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

  if(1 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  PDM_dmesh_t* dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);

  int         *dcell_face_idx;
  PDM_g_num_t *dcell_face;
  int dn_cell = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                           &dcell_face,
                                           &dcell_face_idx,
                                           PDM_OWNERSHIP_KEEP);

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

  PDM_g_num_t* distrib_cell = NULL;
  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_CELL  , &distrib_cell);
  assert(distrib_cell != NULL);

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
  int    *dedge_tag    = (int    * ) malloc(    dn_edge * sizeof(int   ));
  double *dedge_center = (double * ) malloc(3 * dn_edge * sizeof(double));

  for(int i = 0; i < dn_edge; ++i) {

    int i_vtx1 = pedge_vtx[2*i  ]-1;
    int i_vtx2 = pedge_vtx[2*i+1]-1;

    dedge_tag[i] = 0;

    double x1 = pvtx_coord[3*i_vtx1  ];
    double y1 = pvtx_coord[3*i_vtx1+1];
    double z1 = pvtx_coord[3*i_vtx1+2];

    double x2 = pvtx_coord[3*i_vtx2  ];
    double y2 = pvtx_coord[3*i_vtx2+1];
    double z2 = pvtx_coord[3*i_vtx2+2];

    double val1 = _unit_sphere(x1, y1, z1);
    double val2 = _unit_sphere(x2, y2, z2);

    int sgn1 = PDM_SIGN(val1);
    int sgn2 = PDM_SIGN(val2);

    if(sgn1 * sgn2 < 0) {
      dedge_tag[i] = 1;
    }

    dedge_center[3*i  ] = 0.5 * (x1 + x2);
    dedge_center[3*i+1] = 0.5 * (y1 + y2);
    dedge_center[3*i+2] = 0.5 * (z1 + z2);

  }

  int* val = malloc(pn_vtx * sizeof(int));
  printf("pn_vtx = %i \n", pn_vtx);
  for(int i = 0; i < pn_vtx; ++i) {
    double x1 = pvtx_coord[3*i  ];
    double y1 = pvtx_coord[3*i+1];
    double z1 = pvtx_coord[3*i+2];
    val[i]    = PDM_SIGN( _unit_sphere(x1, y1, z1) );
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

  // Compute dcell_edge to get all cell assiciate to an edge that cross iso surf
  int         *dcell_edge_idx = NULL;
  PDM_g_num_t *dcell_edge     = NULL;
  PDM_deduce_combine_connectivity(comm,
                                  distrib_cell,
                                  distrib_face,
                                  dcell_face_idx,
                                  dcell_face,
                                  dface_edge_idx,
                                  dface_edge,
                                  1,
                                  &dcell_edge_idx,
                                  &dcell_edge);

  /*
   * block_to_part on dface_edge
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_edge,
                               (const PDM_g_num_t **) &dcell_edge,
                                                      &dcell_edge_idx[dn_cell],
                                                      1,
                                                      comm);

  int strid_one = 1;
  int **tmp_dcell_edge_tag = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST_INTERLACED,
                         &strid_one,
            (void *  )   dedge_tag,
            (int  ***)   NULL,
            (void ***)  &tmp_dcell_edge_tag);
  int *dcell_edge_tag = tmp_dcell_edge_tag[0];
  free(tmp_dcell_edge_tag);
  free(dedge_tag);

  double **tmp_dcell_edge_center = NULL;
  PDM_block_to_part_exch(btp,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &strid_one,
            (void *  )   dedge_center,
            (int  ***)   NULL,
            (void ***)  &tmp_dcell_edge_center);
  double *dcell_edge_center = tmp_dcell_edge_center[0];
  free(tmp_dcell_edge_center);
  free(dedge_center);

  int         *dcell_tag            = malloc(     dn_cell * sizeof(int        ));
  PDM_g_num_t *cell_to_extract_gnum = malloc(     dn_cell * sizeof(PDM_g_num_t));
  double      *dcell_center         = malloc( 3 * dn_cell * sizeof(double     ));
  int  n_cell_tag = 0;
  int idx_write   = 0;
  for(int i = 0; i < dn_cell; ++i) {
    dcell_tag[i] = 0;



    for(int idx_cell = dcell_edge_idx[i]; idx_cell < dcell_edge_idx[i+1]; ++idx_cell) {
      if(dcell_edge_tag[idx_cell] == 1) {
        dcell_tag[i] = 1;
        cell_to_extract_gnum[n_cell_tag++] = distrib_cell[i_rank] + i + 1;
        break;
      }
    }

    if(dcell_tag[i] == 1) {
      dcell_center[3*idx_write  ] = 0.;
      dcell_center[3*idx_write+1] = 0.;
      dcell_center[3*idx_write+2] = 0.;

      double inv = 1./((double) (dcell_edge_idx[i+1] - dcell_edge_idx[i]));
      for(int idx_cell = dcell_edge_idx[i]; idx_cell < dcell_edge_idx[i+1]; ++idx_cell) {
        dcell_center[3*idx_write  ] += dcell_edge_center[3*idx_cell  ];
        dcell_center[3*idx_write+1] += dcell_edge_center[3*idx_cell+1];
        dcell_center[3*idx_write+2] += dcell_edge_center[3*idx_cell+2];
      }
      dcell_center[3*idx_write  ] = dcell_center[3*idx_write  ] * inv;
      dcell_center[3*idx_write+1] = dcell_center[3*idx_write+1] * inv;
      dcell_center[3*idx_write+2] = dcell_center[3*idx_write+2] * inv;

      idx_write++;
    }

  }
  free(dcell_edge_center);

  if(0 == 1) {
    PDM_log_trace_array_int (dedge_tag, dn_cell, "dcell_tag");
    PDM_log_trace_array_long(cell_to_extract_gnum, n_cell_tag, "cell_to_extract_gnum");
  }

  PDM_block_to_part_free(btp);

  free(edge_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(pedge_vtx_idx);
  free(pedge_vtx);
  free(dcell_edge_tag);

  // Extract

  sprintf(filename, "out_equi_cell_coord_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            n_cell_tag,
                            dcell_center,
                            NULL,
                            NULL);



  /*
   * Rebuild partition that contains faces and reequilibrate
   */
  PDM_gen_gnum_t* gnum_equi = PDM_gnum_create(3, 1, PDM_FALSE, 0., comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_coords(gnum_equi, 0, n_cell_tag, dcell_center, NULL);
  PDM_gnum_compute(gnum_equi);
  PDM_g_num_t* child_equi_cell_gnum = PDM_gnum_get(gnum_equi, 0);
  PDM_gnum_free(gnum_equi);
  free(dcell_center);

  /*
   * Equilibrage avec le part_to_block
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                       1.,
                                                       &child_equi_cell_gnum,
                                                       NULL,
                                                       &n_cell_tag,
                                                       1,
                                                       comm);

  int n_cell_equi = PDM_part_to_block_n_elt_block_get (ptb);
  PDM_g_num_t *block_cell_g_num_child_equi = PDM_part_to_block_block_gnum_get (ptb);

  PDM_g_num_t *block_cell_equi_parent_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
               (void **) &cell_to_extract_gnum,
                          NULL,
               (void **) &block_cell_equi_parent_g_num);

  if(0 == 1) {
    PDM_log_trace_array_long(block_cell_equi_parent_g_num, n_cell_equi, "block_cell_equi_parent_g_num ::");
    PDM_log_trace_array_long(block_cell_g_num_child_equi , n_cell_equi, "block_cell_g_num_child_equi  ::");
  }

  /*
   * Je prepare tout pour mon petit Bastien
   */
  int          pn_edge_equi        = 0;
  PDM_g_num_t *pequi_parent_edge_ln_to_gn = NULL;
  int         *pequi_cell_edge_idx = NULL;
  int         *pequi_cell_edge     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_cell,
                                                           dcell_edge_idx,
                                                           dcell_edge,
                                                           n_cell_equi,
                                     (const PDM_g_num_t *) block_cell_equi_parent_g_num,
                                                           &pn_edge_equi,
                                                           &pequi_parent_edge_ln_to_gn,
                                                           &pequi_cell_edge_idx,
                                                           &pequi_cell_edge);

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

  _iso_surf(comm,
            n_cell_equi,
            pn_edge_equi,
            pn_vtx_equi,
            pequi_cell_edge_idx,
            pequi_cell_edge,
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
  free(pequi_cell_edge_idx);
  free(pequi_cell_edge);


  free(block_cell_equi_parent_g_num);


  PDM_part_to_block_free(ptb);
  free(dcell_edge_idx);
  free(dcell_edge);
  free(child_equi_cell_gnum);


  free(cell_to_extract_gnum);
  free(dcell_tag);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_dcube_nodal_gen_free(dcube);


  PDM_MPI_Finalize();

  return 0;
}
