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
#include "pdm_poly_vol_gen.h"
#include "pdm_distrib.h"
#include "pdm_multipart.h"
#include "pdm_iso_surface.h"
#include "pdm_iso_surface_priv.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_priv.h"
#include "pdm_writer.h"



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





#define PERLIN_NOISE_N 8
static double perlin_fx[PERLIN_NOISE_N][PERLIN_NOISE_N][PERLIN_NOISE_N];
static double perlin_fy[PERLIN_NOISE_N][PERLIN_NOISE_N][PERLIN_NOISE_N];
static double perlin_fz[PERLIN_NOISE_N][PERLIN_NOISE_N][PERLIN_NOISE_N];
static const double perlin_step = 1. / (double) (PERLIN_NOISE_N - 1);

static void
_init_perlin_noise (void)
{
  double i_rand_max = 1./(double) RAND_MAX;

  for (int k = 0; k < PERLIN_NOISE_N; k++) {
    for (int j = 0; j < PERLIN_NOISE_N; j++) {
      for (int i = 0; i < PERLIN_NOISE_N; i++) {
        perlin_fx[i][j][k] = 2*rand()*i_rand_max - 1;
        perlin_fy[i][j][k] = 2*rand()*i_rand_max - 1;
        perlin_fz[i][j][k] = 2*rand()*i_rand_max - 1;
      }
    }
  }
}

static void
_update_perlin_noise (const double dt)
{
  const double c = cos(dt);
  const double s = sin(dt);

  for (int k = 0; k < PERLIN_NOISE_N; k++) {
    for (int j = 0; j < PERLIN_NOISE_N; j++) {
      for (int i = 0; i < PERLIN_NOISE_N; i++) {
        double x = perlin_fx[i][j][k];
        double y = perlin_fy[i][j][k];
        double z = perlin_fz[i][j][k];

        perlin_fx[i][j][k] =  c*x + s*y;
        perlin_fy[i][j][k] = -s*x + c*y;
        perlin_fz[i][j][k] = z;
      }
    }
  }
}


static double _smooth_step (const double x)
{
  return x*(3*x - 2*x*x);
}

static double _lerp (const double a, const double b, const double x)
{
  return (1 - x)*a + x*b;
}

static double _dot (
                    const double x, const double y, const double z,
                    const int i, const int j, const int k
                    )
{
  return
  x*perlin_fx[i][j][k] +
  y*perlin_fy[i][j][k] +
  z*perlin_fz[i][j][k];
}

static double
_eval_perlin
(
 const double  x,
 const double  y,
 const double  z
 )
 {
  // #define _smooth_step (x) ((x)*(3*(x) - 2*(x)*(x)))
  // #define _lerp (a,b,x) ((1 - (x))*(a) + (x)*(b))
  // #define _dot (x,y,z,i,j,k) ((x)*perlin_fx[(i)][(j)][(k)] + (y)*perlin_fy[(i)][(j)][(k)] + (z)*perlin_fz[(i)][(j)][(k)])

  int i = (int) (x / perlin_step);
  int j = (int) (y / perlin_step);
  int k = (int) (z / perlin_step);

  i = PDM_MIN(PDM_MAX(i, 0), PERLIN_NOISE_N-2);
  j = PDM_MIN(PDM_MAX(j, 0), PERLIN_NOISE_N-2);
  k = PDM_MIN(PDM_MAX(k, 0), PERLIN_NOISE_N-2);

  double tx = x - i*perlin_step;
  double ty = y - j*perlin_step;
  double tz = z - k*perlin_step;

  double u = _smooth_step(tx/perlin_step);
  double v = _smooth_step(ty/perlin_step);
  double w = _smooth_step(tz/perlin_step);

  double x0 = tx;
  double y0 = ty;
  double z0 = tz;
  double x1 = tx - perlin_step;
  double y1 = ty - perlin_step;
  double z1 = tz - perlin_step;


  double d000 = _dot(x0, y0, z0, i,   j,   k  );
  double d100 = _dot(x1, y0, z0, i+1, j,   k  );
  double d010 = _dot(x0, y1, z0, i,   j+1, k  );
  double d110 = _dot(x1, y1, z0, i+1, j+1, k  );
  double d001 = _dot(x0, y0, z1, i,   j,   k+1);
  double d101 = _dot(x1, y0, z1, i+1, j,   k+1);
  double d011 = _dot(x0, y1, z1, i,   j+1, k+1);
  double d111 = _dot(x1, y1, z1, i+1, j+1, k+1);

  double a = _lerp(d000, d100, u);
  double b = _lerp(d010, d110, u);
  double c = _lerp(d001, d101, u);
  double d = _lerp(d011, d111, u);

  double e = _lerp(a, b, v);
  double f = _lerp(c, d, v);

  double val = _lerp(e, f, w) - 0.041;

  // #undef _smooth_step
  // #undef _lerp
  // #undef _dot

  return val;
 }

static void
_eval_perlin_noise
(
 const double  x,
 const double  y,
 const double  z,
       double *f,
       double *df_dx,
       double *df_dy,
       double *df_dz
 )
{
  const double h = 1e-4;

  *f = _eval_perlin(x, y, z);
  *df_dx = (_eval_perlin(x + h, y,     z    ) - (*f)) / h;
  *df_dy = (_eval_perlin(x,     y + h, z    ) - (*f)) / h;
  *df_dz = (_eval_perlin(x,     y,     z + h) - (*f)) / h;
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
           int                   *n_part,
           int                   *part_method,
           int                   *n_steps,
           double                *dt)
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
    else if (strcmp(argv[i], "-n_steps") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *n_steps = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-dt") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *dt = atof(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}







static void
_write_instant
(
 const int          step,
 PDM_iso_surface_t *isos
 )
{
  int     n_vtx;
  int     n_face;
  int    *face_vtx;
  double *vtx_coord;
  int    *face_rank;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);
  PDM_MPI_Comm_size(isos->comm, &n_rank);


  if (n_rank == 1) {
    n_face    = isos->isosurf_n_face;
    n_vtx     = isos->isosurf_n_vtx;
    face_vtx  = isos->isosurf_face_vtx;
    vtx_coord = isos->isosurf_vtx_coord;

    face_rank = PDM_array_zeros_int(n_face);
  }

  else {

    PDM_g_num_t lmax_face_gnum = 0;
    for (int i = 0; i < isos->isosurf_n_face; i++) {
      lmax_face_gnum = PDM_MAX(lmax_face_gnum, isos->isosurf_face_ln_to_gn[i]);
    }

    PDM_g_num_t lmax_vtx_gnum = 0;
    for (int i = 0; i < isos->isosurf_n_vtx; i++) {
      lmax_vtx_gnum = PDM_MAX(lmax_vtx_gnum, isos->isosurf_vtx_ln_to_gn[i]);
    }

    PDM_g_num_t gn_face;
    PDM_MPI_Allreduce(&lmax_face_gnum, &gn_face, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, isos->comm);

    PDM_g_num_t gn_vtx;
    PDM_MPI_Allreduce(&lmax_vtx_gnum ,  &gn_vtx,  1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, isos->comm);

    PDM_g_num_t *distrib_face = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_rank + 1));
    PDM_g_num_t *distrib_vtx  = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_rank + 1));
    distrib_face[0] = 0;
    distrib_vtx [0] = 0;
    for (int i = 1; i <= n_rank; i++) {
      distrib_face[i] = gn_face;
      distrib_vtx [i] = gn_vtx;
    }

    // PDM_log_trace_array_long(distrib_face, n_rank+1, "distrib_face : ");
    // PDM_log_trace_array_long(distrib_vtx,  n_rank+1, "distrib_vtx : ");



    int *pface_vtx = (int *) malloc(sizeof(int) * 3*isos->isosurf_n_face);
    for (int i = 0; i < 3*isos->isosurf_n_face; i++) {
      pface_vtx[i] = (int) isos->isosurf_vtx_ln_to_gn[isos->isosurf_face_vtx[i]-1];
    }


    PDM_part_to_block_t *ptb_face = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                          PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                          1.,
                                                                          &isos->isosurf_face_ln_to_gn,
                                                    (const PDM_g_num_t *) distrib_face,
                                                                          &isos->isosurf_n_face,
                                                                          1,
                                                                          isos->comm);

    n_face = PDM_part_to_block_n_elt_block_get(ptb_face);

    PDM_part_to_block_exch(ptb_face,
                           sizeof(int),
                           PDM_STRIDE_CST_INTERLACED,
                           3,
                           NULL,
                 (void **) &pface_vtx,
                           NULL,
                 (void **) &face_vtx);
    free(pface_vtx);


    int *pface_rank = (int *) malloc(sizeof(int) * isos->isosurf_n_face);
    for (int i = 0; i < isos->isosurf_n_face; i++) {
      pface_rank[i] = (int) i_rank;
    }

    PDM_part_to_block_exch(ptb_face,
                           sizeof(int),
                           PDM_STRIDE_CST_INTERLACED,
                           1,
                           NULL,
                 (void **) &pface_rank,
                           NULL,
                 (void **) &face_rank);
    free(pface_rank);

    PDM_part_to_block_free(ptb_face);




    PDM_part_to_block_t *ptb_vtx = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                         PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                         1.,
                                                                         &isos->isosurf_vtx_ln_to_gn,
                                                   (const PDM_g_num_t *) distrib_vtx,
                                                                         &isos->isosurf_n_vtx,
                                                                         1,
                                                                         isos->comm);

    n_vtx = PDM_part_to_block_n_elt_block_get(ptb_vtx);

    PDM_part_to_block_exch(ptb_vtx,
                           3*sizeof(double),
                           PDM_STRIDE_CST_INTERLACED,
                           1,
                           NULL,
                 (void **) &isos->isosurf_vtx_coord,
                           NULL,
                 (void **) &vtx_coord);

    PDM_part_to_block_free(ptb_vtx);

    free(distrib_face);
    free(distrib_vtx);

  }

  if (i_rank == 0) {
    const char* field_name[] = {"num_part", 0 };

    char filename[999];
    sprintf(filename, "isosurf_unsteady_%4.4d.vtk", step);
    PDM_vtk_write_std_elements(filename,
                               n_vtx,
                               vtx_coord,
                               NULL,
                               PDM_MESH_NODAL_TRIA3,
                               n_face,
                               face_vtx,
                               NULL,
                               1,
                               field_name,
                (const int **) &face_rank);
    free(face_rank);
  }


  if (n_rank > 1 && i_rank == 1) {
    free(face_vtx);
    free(vtx_coord);
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

  PDM_g_num_t          n_vtx_seg = 10;
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_TETRA4;
  //  5 -> tetra
  //  6 -> pyramid
  //  7 -> prism
  //  8 -> hexa
  //  9 -> poly3d
  int                  n_part    = -1;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  int    n_steps = 2;
  double dt      = 0.1;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &elt_type,
             &n_part,
     (int *) &part_method,
             &n_steps,
             &dt);
  assert(PDM_Mesh_nodal_elt_dim_get(elt_type) == 3);
  assert(n_part < 1);

  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_multipart_t *mpart = NULL;
  int n_zone = 1;
  int *n_part_zones = NULL;
  if (n_part > 0) {
    n_part_zones = (int *) malloc(sizeof(int) * n_zone);
    n_part_zones[0] = n_part;
    mpart = PDM_multipart_create(n_zone,
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
  }


  /*
   *  Create distributed cube
   */
  PDM_dmesh_t *dmesh = NULL;

  PDM_g_num_t *distrib_cell = NULL;
  PDM_g_num_t *distrib_face = NULL;
  PDM_g_num_t *distrib_edge = NULL;
  PDM_g_num_t *distrib_vtx  = NULL;

  int          dn_cell = 0;
  int          dn_face = 0;
  int          dn_edge = 0;
  int          dn_vtx  = 0;
  double      *dvtx_coord      = NULL;
  int         *dcell_face_idx  = NULL;
  PDM_g_num_t *dcell_face      = NULL;
  int         *dface_edge_idx  = NULL;
  PDM_g_num_t *dface_edge      = NULL;
  int         *dedge_vtx_idx   = NULL;
  PDM_g_num_t *dedge_vtx       = NULL;

  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         1.,
                                                         0.,
                                                         0.,
                                                         0.,
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];

  PDM_dmesh_nodal_to_dmesh_t *dmntodm = PDM_dmesh_nodal_to_dmesh_create(1,
                                                                        comm,
                                                                        PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);

  dn_cell = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                       &dcell_face,
                                       &dcell_face_idx,
                                       PDM_OWNERSHIP_KEEP);

  dn_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                       &dface_edge,
                                       &dface_edge_idx,
                                       PDM_OWNERSHIP_KEEP);

  dn_edge = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                       &dedge_vtx,
                                       &dedge_vtx_idx,
                                       PDM_OWNERSHIP_KEEP);

  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_EDGE, &distrib_edge);
  assert(distrib_edge != NULL);

  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE, &distrib_face);
  assert(distrib_face != NULL);

  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_CELL, &distrib_cell);
  assert(distrib_cell != NULL);


  if (n_part > 0) {
    PDM_multipart_register_dmesh_nodal(mpart, 0, dmn);

    PDM_multipart_run_ppart(mpart);
  }




  _init_perlin_noise();

  PDM_iso_surface_t* isos = PDM_iso_surface_create(3,
                                                   PDM_ISO_SURFACE_KIND_FIELD,
                                                   1,
                                                   PDM_OWNERSHIP_KEEP,
                                                   comm);




  /*
   *  Allocate field and gradient and set isosurface input
   */

  double  *dfield          = NULL;
  double  *dgradient_field = NULL;
  double **pfield          = NULL;
  double **pgradient_field = NULL;

  if (n_part > 0) {
    pfield          = (double **) malloc(sizeof(double *) * n_part);
    pgradient_field = (double **) malloc(sizeof(double *) * n_part);

    for (int i_part = 0; i_part < n_part; i_part++) {
      double *vtx_coord;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx_coord,
                                                   PDM_OWNERSHIP_KEEP);
      pfield[i_part]          = (double *) malloc(sizeof(double) * n_vtx);
      pgradient_field[i_part] = (double *) malloc(sizeof(double) * n_vtx * 3);
    }

    // PDM_iso_surface_part_set (isos,
    //                             i_part,
    //                             n_cell,
    //                             n_face,
    //                             n_edge,
    //                             n_vtx,
    //                             cell_face_idx,
    //                             cell_face,
    //                             face_edge_idx,
    //                             face_edge,
    //                             edge_vtx,
    //                             NULL,//face_vtx_idx,
    //                             NULL,//face_vtx,
    //                             cell_ln_to_gn,
    //                             face_ln_to_gn,
    //                             edge_ln_to_gn,
    //                             vtx_ln_to_gn,
    //                             vtx_coord);

    //   PDM_iso_surface_part_field_set (isos,
    //                                   i_part,
    //                                   pfield[i_part]);

    //   PDM_iso_surface_part_gradient_field_set (isos,
    //                                            i_part,
    //                                            pgradient_field[i_part]);
  }

  else {

    dfield          = (double *) malloc(    dn_vtx * sizeof(double));
    dgradient_field = (double *) malloc(3 * dn_vtx * sizeof(double));

    PDM_iso_surface_dconnectivity_set(isos,
                                      PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                      dcell_face,
                                      dcell_face_idx);
    PDM_iso_surface_dconnectivity_set(isos,
                                      PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                      dface_edge,
                                      dface_edge_idx);
    PDM_iso_surface_dconnectivity_set(isos,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      dedge_vtx,
                                      dedge_vtx_idx);

    PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_CELL  , distrib_cell);
    PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_FACE  , distrib_face);
    PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_EDGE  , distrib_edge);
    PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_VERTEX, distrib_vtx);


    PDM_iso_surface_dvtx_coord_set (isos, dvtx_coord     );
    PDM_iso_surface_dfield_set     (isos, dfield         );
    // PDM_iso_surface_dgrad_field_set(isos, dgradient_field);

  }




  /*
   *  Main loop
   */

  // PDM_writer_t *cs = PDM_writer_create ("Ensight",
  //                                       PDM_WRITER_FMT_ASCII,
  //                                       PDM_WRITER_TOPO_VARIABLE,
  //                                       PDM_WRITER_OFF,
  //                                       "isosurf_unsteady",
  //                                       "isosurf_unsteady",
  //                                       comm,
  //                                       PDM_IO_ACCES_MPI_SIMPLE,
  //                                       1.,
  //                                       NULL);

  // int id_geom = PDM_writer_geom_create (cs,
  //                                       "isosurf_unsteady",
  //                                       PDM_WRITER_OFF,
  //                                       PDM_WRITER_OFF,
  //                                       1);

  // int id_var_part = PDM_writer_var_create (cs,
  //                                          PDM_WRITER_OFF,
  //                                          PDM_WRITER_VAR_SCALAIRE,
  //                                          PDM_WRITER_VAR_ELEMENTS,
  //                                          "num_part");

  // int id_block = PDM_writer_geom_bloc_add(cs,
  //                                         id_geom,
  //                                         PDM_WRITER_ON,
  //                                         PDM_WRITER_TRIA3);

  for (int step = 0; step < n_steps; step++) {

    if (i_rank == 0) {
      printf("\nStep %d/%d\n", step+1, n_steps);
      fflush(stdout);
    }

    /*
     *  Reset iso-surface
     */
    isos->isosurf_n_vtx  = 0;
    isos->isosurf_n_face = 0;
    if (isos->isosurf_face_vtx_idx  != NULL) free(isos->isosurf_face_vtx_idx );
    if (isos->isosurf_face_vtx      != NULL) free(isos->isosurf_face_vtx     );
    if (isos->isosurf_vtx_coord     != NULL) free(isos->isosurf_vtx_coord    );
    if (isos->isosurf_face_ln_to_gn != NULL) free(isos->isosurf_face_ln_to_gn);
    if (isos->isosurf_vtx_ln_to_gn  != NULL) free(isos->isosurf_vtx_ln_to_gn );


    /*
     *  Update field and gradient
     */
     _update_perlin_noise(dt);

    if (n_part > 0) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        double *vtx_coord;
        int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                     0,
                                                     i_part,
                                                     &vtx_coord,
                                                     PDM_OWNERSHIP_KEEP);
        for (int i = 0; i < n_vtx; i++) {
          _eval_perlin_noise(vtx_coord[3*i  ],
                             vtx_coord[3*i+1],
                             vtx_coord[3*i+2],
                             &pfield[i_part][i],
                             &pgradient_field[i_part][3*i],
                             &pgradient_field[i_part][3*i+1],
                             &pgradient_field[i_part][3*i+2]);
        }
      }
    }

    else {

     for (int i = 0; i < dn_vtx; i++) {
       _eval_perlin_noise(dvtx_coord[3*i  ],
                          dvtx_coord[3*i+1],
                          dvtx_coord[3*i+2],
                          &dfield[i],
                          &dgradient_field[3*i],
                          &dgradient_field[3*i+1],
                          &dgradient_field[3*i+2]);
     }
   }


    /*
     *  Compute iso-surface
     */
    PDM_iso_surface_compute(isos);


    /*
     *  Write new step
     */
    _write_instant(step, isos);
    // int          isosurf_n_vtx;
    // int          isosurf_n_face;
    // int         *isosurf_face_vtx_idx  = NULL;
    // int         *isosurf_face_vtx      = NULL;
    // double      *isosurf_vtx_coord     = NULL;
    // PDM_g_num_t *isosurf_face_ln_to_gn = NULL;
    // PDM_g_num_t *isosurf_vtx_ln_to_gn  = NULL;
    // PDM_iso_surface_surface_get (isos,
    //                              &isosurf_n_vtx,
    //                              &isosurf_n_face,
    //                              &isosurf_face_vtx_idx,
    //                              &isosurf_face_vtx,
    //                              &isosurf_vtx_coord,
    //                              &isosurf_face_ln_to_gn,
    //                              &isosurf_vtx_ln_to_gn);
    // printf("vtx1: %p %p\n", (void *) isos->isosurf_vtx_coord, (void *) isos->isosurf_vtx_ln_to_gn);

    // PDM_writer_step_beg (cs, step*dt);


    // /* Write geometry */
    // PDM_writer_geom_bloc_std_set(cs,
    //                              id_geom,
    //                              id_block,
    //                              0,
    //                              isos->isosurf_n_face,
    //                              isos->isosurf_face_vtx,
    //                              isos->isosurf_face_ln_to_gn);
    // PDM_writer_geom_coord_set(cs,
    //                           id_geom,
    //                           0,
    //                           isos->isosurf_n_vtx,
    //                           isos->isosurf_vtx_coord,
    //                           isos->isosurf_vtx_ln_to_gn);
    // log_trace("step %d\n", step+1);
    // PDM_log_trace_array_long(isos->isosurf_vtx_ln_to_gn,
    //                          isos->isosurf_n_vtx,
    //                          "isosurf_vtx_ln_to_gn : ");
    // PDM_log_trace_connectivity_int(isos->isosurf_face_vtx_idx,
    //                                isos->isosurf_face_vtx,
    //                                isos->isosurf_n_face,
    //                                "isosurf_face_vtx : ");


    // PDM_writer_geom_write(cs,
    //                       id_geom);

    // /* Write variables */
    // PDM_real_t *val_part = (PDM_real_t *) malloc(sizeof(PDM_real_t) * isos->isosurf_n_face);
    // for (int i = 0; i < isos->isosurf_n_face; i++) {
    //   val_part[i] = (PDM_real_t) i_rank;
    // }

    // PDM_writer_var_set(cs,
    //                    id_var_part,
    //                    id_geom,
    //                    0,
    //                    (const PDM_real_t *) val_part);

    // PDM_writer_var_write(cs,
    //                      id_var_part);

    // PDM_writer_step_end(cs);
    // free(val_part);

    // PDM_writer_geom_data_free(cs, id_geom);
    // printf("vtx2: %p %p\n", (void *) isos->isosurf_vtx_coord, (void *) isos->isosurf_vtx_ln_to_gn);
    // isos->isosurf_face_vtx_idx  = NULL;
    // isos->isosurf_face_vtx      = NULL;
    // isos->isosurf_vtx_coord     = NULL;
    // isos->isosurf_face_ln_to_gn = NULL;
    // isos->isosurf_vtx_ln_to_gn  = NULL;

  } // main loop

  // PDM_writer_free(cs);

  /*
   *  Free memory
   */
  if (n_part > 0) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      free(pfield[i_part]);
      free(pgradient_field[i_part]);
    }
    free(pfield);
    free(pgradient_field);

    PDM_multipart_free(mpart);
  }
  else {
    free(dfield);
    free(dgradient_field);
  }

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_dcube_nodal_gen_free(dcube);

  isos->isosurf_face_vtx      = NULL;
  isos->isosurf_face_ln_to_gn = NULL;
  PDM_iso_surface_free(isos);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();
}
