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
#include "pdm_poly_surf_gen.h"
#include "pdm_distrib.h"
#include "pdm_iso_surface.h"

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
           double                *length,
           PDM_Mesh_nodal_elt_t  *elt_type)
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
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static const double _A = 1.072;
static const double _B = 1.044;
static const double _C = 0.286;
static const double _D = -0.042;
static const double _E = 0.067;
static const double _F = -0.025;


static
inline
double
_unit_circle
(
 double x,
 double y
)
{
  // return x * x + y * y - 0.1;
  // return _A*x*x + _B*x*y + _C*y*y + _D*x + _E*y + _F;
  // double a = (x*x + y*y - 1);
  // return a*a*a - x*x*(6*y*x*x + x*x - 2*y*y*y - 6*y);
  // return a*a*a - x*x*(x*x*(1 + 6*y) - y*(6 + 2*y*y));


  // return x*x*x*x*x*x - y*x*(x*x-y*y) + y*y*y*y*y*y;
  // return PDM_MAX(PDM_ABS(x), PDM_ABS(y)) - 0.25;
  // return PDM_MIN(PDM_ABS(x) + PDM_ABS(y) - 0.25, (x-0.1)*(x-0.1) + (y+0.2)*(y+0.2) - 0.07);

  double v1 = (x-0.23)*(x-0.23) + (y-0.28)*(y-0.28) - 0.03;
  double v2 = (x+0.23)*(x+0.23) + (y-0.28)*(y-0.28) - 0.03;
  double v3 = x*x + y*y - 0.1;
  return PDM_MIN(PDM_MIN(v1, v2), v3);
}


static
inline
void
_unit_circle_gradient
(
 double  x,
 double  y,
 double *df_dx,
 double *df_dy
)
{
  *df_dx = 2*x;
  *df_dy = 2*y;
  // return;
  *df_dx = 2*_A*x + _B*y + _D;
  *df_dy = 2*_C*y + _B*x + _E;

  *df_dx = 6*x*x*x*x*x + (12*y*y - 24*y - 16)*x*x*x + (6*y*y*y*y + 4*y*y*y - 12*y*y + 12*y + 6)*x;
  // *df_dx = x*( (6 + 2*y*(6 + y*(6 + y*(2 + 3*y)))) + x*x*(16 + 12*y*(-2 + y)) + 6*x*x );
  double a = (x*x + y*y - 1);
  *df_dy = 6*y*a*a - x*x*(-6*y*y + 6*x*x - 6);
  // *df_dy = 6*(y*a*a - x*x*(-y*y + x*x - 1));


  *df_dx = 6*x*x*x*x*x - 3*y*x*x + y*y*y;
  *df_dx = 6*y*y*y*y*y + 3*x*y*y - x*x*x;


  // if (PDM_ABS(x) > PDM_ABS(y)) {
  //   *df_dx = PDM_SIGN(x);
  //   *df_dy = 0.;
  // } else {
  //   *df_dx = 0;
  //   *df_dy = PDM_SIGN(y);
  // }
  // if (PDM_ABS(x) + PDM_ABS(y) - 0.25 < (x-0.1)*(x-0.1) + (y+0.2)*(y+0.2) - 0.07) {
  //   *df_dx = PDM_SIGN(x);
  //   *df_dy = PDM_SIGN(y);
  // } else {
  //   *df_dx = 2*(x-0.1);
  //   *df_dy = 2*(y+0.2);
  // }

  double df_dx1 = 2*(x-0.23);
  double df_dy1 = 2*(y-0.28);
  double df_dx2 = 2*(x+0.23);
  double df_dy2 = 2*(y-0.28);
  double df_dx3 = 2*x;
  double df_dy3 = 2*y;

  double v1 = (x-0.23)*(x-0.23) + (y-0.28)*(y-0.28) - 0.03;
  double v2 = (x+0.23)*(x+0.23) + (y-0.28)*(y-0.28) - 0.03;
  double v3 = x*x + y*y - 0.1;

  if (v1 < v2) {
    if (v1 < v3) {
      *df_dx = df_dx1;
      *df_dy = df_dy1;
    } else {
      *df_dx = df_dx3;
      *df_dy = df_dy3;
    }
  } else {
    if (v2 < v3) {
      *df_dx = df_dx2;
      *df_dy = df_dy2;
    } else {
      *df_dx = df_dx3;
      *df_dy = df_dy3;
    }
  }
}



static
void
_dump_dmesh
(
 PDM_MPI_Comm  comm,
 const int     dn_face,
 const int     dn_vtx,
 int          *dface_vtx_idx,
 PDM_g_num_t  *dface_vtx,
 double       *dvtx_coord
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t *distrib_face = PDM_compute_entity_distribution(comm, dn_face);
  PDM_g_num_t *distrib_vtx  = PDM_compute_entity_distribution(comm, dn_vtx);

  PDM_g_num_t *dface_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dn_face);
  for (int i = 0; i < dn_face; i++) {
    dface_ln_to_gn[i] = distrib_face[i_rank] + i + 1;
  }

  int          pn_vtx        = 0;
  PDM_g_num_t *pvtx_ln_to_gn = NULL;
  int         *pface_vtx_idx = NULL;
  int         *pface_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_vtx_idx,
                                                           dface_vtx,
                                                           dn_face,
                                     (const PDM_g_num_t *) dface_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pface_vtx_idx,
                                                           &pface_vtx);

  double** tmp_pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        distrib_vtx,
                                        dvtx_coord,
                                        &pn_vtx,
                 (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coord);
  double* pvtx_coord = tmp_pvtx_coord[0];
  free (tmp_pvtx_coord);

  char filename[999];
  sprintf(filename, "mesh_%2.2d.vtk", i_rank);
  PDM_vtk_write_polydata(filename,
                         pn_vtx,
                         pvtx_coord,
                         pvtx_ln_to_gn,
                         dn_face,
                         pface_vtx_idx,
                         pface_vtx,
                         dface_ln_to_gn,
                         NULL);

  free(pvtx_ln_to_gn);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pvtx_coord);
  free(dface_ln_to_gn);
  free(distrib_face);
  free(distrib_vtx);
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
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_TRIA3;
  //  2 -> tria
  //  3 -> quad
  //  4 -> poly2d

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &elt_type);

  assert(PDM_Mesh_nodal_elt_dim_get(elt_type) == 2);

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
   *  Create distributed cube
   */

  PDM_g_num_t *distrib_face = NULL;
  PDM_g_num_t *distrib_edge = NULL;
  PDM_g_num_t *vtx_distrib  = NULL;

  PDM_dcube_nodal_t          *dcube = NULL;
  PDM_dmesh_nodal_to_dmesh_t *dmntodm = NULL;

  PDM_g_num_t     ng_face         = 0;
  PDM_g_num_t     ng_vtx          = 0;
  PDM_g_num_t     ng_edge         = 0;
  int             dn_vtx          = 0;
  int             dn_face         = 0;
  int             dn_edge         = 0;
  int             n_edge_group    = 0;
  double         *dvtx_coord      = NULL;
  int            *dface_edge_idx  = NULL;
  PDM_g_num_t    *dface_vtx       = NULL;
  PDM_g_num_t    *dface_edge      = NULL;
  int            *dedge_vtx_idx   = NULL;
  PDM_g_num_t    *dedge_vtx       = NULL;
  PDM_g_num_t    *dedge_face      = NULL;
  int            *dedge_group_idx = NULL;
  PDM_g_num_t    *dedge_group     = NULL;


  if (elt_type < PDM_MESH_NODAL_POLY_2D) {

    dcube = PDM_dcube_nodal_gen_create (comm,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        length,
                                        -0.5*length,// -0.35,
                                        -0.5*length,// -0.3,
                                        0.,
                                        elt_type,
                                        1,
                                        PDM_OWNERSHIP_KEEP);
    PDM_dcube_nodal_gen_build (dcube);


    PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
    PDM_dmesh_nodal_generate_distribution(dmn);

    vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];


    if (0) {
      double noise = 0.2 * length / (double) (n_vtx_seg - 1);
      for (int i = 0; i < dn_vtx; i++) {
        for (int j = 0; j < 2; j++) {
          dvtx_coord[3*i+j] += noise * (rand() / (double) RAND_MAX - 0.5);
        }
      }
    }

    if(1 == 1) {
    // PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
    }

    dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

    PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

    PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

    PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);

    PDM_dmesh_t* dmesh = NULL;
    PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);


    dn_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                         &dface_edge,
                                         &dface_edge_idx,
                                         PDM_OWNERSHIP_KEEP);

    dn_edge  = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &dedge_vtx,
                                          &dedge_vtx_idx,
                                          PDM_OWNERSHIP_KEEP);

    if(0 == 1) {
      PDM_log_trace_connectivity_long(dface_edge_idx, dface_edge, dn_face, "dface_edge ::");
      PDM_log_trace_connectivity_long(dedge_vtx_idx , dedge_vtx , dn_edge, "dedge_vtx  ::");
    }

    PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_EDGE  , &distrib_edge);
    assert(distrib_edge != NULL);

    PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE  , &distrib_face);
    assert(distrib_face != NULL);
  }

  else {

    PDM_poly_surf_gen (comm,
                       -0.5*length,
                       0.5*length,
                       -0.5*length,
                       0.5*length,
                       1,
                       0,
                       n_vtx_seg,
                       n_vtx_seg,
                       &ng_face,
                       &ng_vtx,
                       &ng_edge,
                       &dn_vtx,
                       &dvtx_coord,
                       &dn_face,
                       &dface_edge_idx,
                       &dface_vtx,
                       &dface_edge,
                       &dn_edge,
                       &dedge_vtx,
                       &dedge_face,
                       &n_edge_group,
                       &dedge_group_idx,
                       &dedge_group);

    _dump_dmesh (comm,
                 dn_face,
                 dn_vtx,
                 dface_edge_idx,
                 dface_vtx,
                 dvtx_coord);

    distrib_face = PDM_compute_entity_distribution(comm, dn_face);
    distrib_edge = PDM_compute_entity_distribution(comm, dn_edge);
    vtx_distrib  = PDM_compute_entity_distribution(comm, dn_vtx);

    dedge_vtx_idx = (int *) malloc(sizeof(int) * (dn_edge + 1));
    for (int i = 0; i <= dn_edge; i++) {
      dedge_vtx_idx[i] = 2*i;
    }
  }

  // Compute dfield and gradient field
  double *dfield          = (double *) malloc(     dn_vtx * sizeof(double));
  double *dgradient_field = (double *) malloc( 3 * dn_vtx * sizeof(double));
  // double *dgradient_field = NULL;

  for(int i = 0; i < dn_vtx; ++i) {

    double x1 = dvtx_coord[3*i  ];
    double y1 = dvtx_coord[3*i+1];
    dfield[i] = _unit_circle(x1, y1);

    _unit_circle_gradient(x1, y1, &dgradient_field[3*i], &dgradient_field[3*i+1]);

    dgradient_field[3*i+2] = 0;
  }

  // Taylor Green
  // for(int i = 0; i < dn_vtx; ++i) {

  //   double x1 = dvtx_coord[3*i  ];
  //   double y1 = dvtx_coord[3*i+1];
  //   dfield[i] = _unit_circle(x1, y1);

  //   _unit_circle_gradient(x1, y1, &dgradient_field[3*i], &dgradient_field[3*i+1]);

  //   dgradient_field[3*i+2] = 0;
  // }

  PDM_iso_surface_t* isos = PDM_iso_surface_create(2, PDM_ISO_SURFACE_KIND_FIELD, 1, PDM_OWNERSHIP_KEEP, comm);
  // PDM_iso_surface_t* isos = PDM_iso_surface_create(2, PDM_ISO_SURFACE_KIND_PLANE, 1, PDM_OWNERSHIP_KEEP, comm);

  PDM_iso_surface_plane_equation_set(isos, 1., 0., 0., -0);
  // PDM_iso_surface_plane_equation_set(isos, 1., 0.5, 0.25, -0.0234);

  PDM_iso_surface_dconnectivity_set(isos,
                                    PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                    dface_edge,
                                    dface_edge_idx);
  PDM_iso_surface_dconnectivity_set(isos,
                                    PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                    dedge_vtx,
                                    dedge_vtx_idx);

  PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_FACE  , distrib_face);
  PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_EDGE  , distrib_edge);
  PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_VERTEX, vtx_distrib);

  PDM_iso_surface_dvtx_coord_set (isos, dvtx_coord     );
  PDM_iso_surface_dfield_set     (isos, dfield         );
  PDM_iso_surface_dgrad_field_set(isos, dgradient_field);

  PDM_iso_surface_compute(isos);

  PDM_iso_surface_free(isos);

  free(dfield);
  free(dgradient_field);

  if (elt_type < PDM_MESH_NODAL_POLY_2D) {
    PDM_dmesh_nodal_to_dmesh_free(dmntodm);
    PDM_dcube_nodal_gen_free(dcube);
  } else {
    free(distrib_face);
    free(distrib_edge);
    free(vtx_distrib);
    free(dvtx_coord);
    free(dface_edge_idx);
    free(dface_vtx);
    free(dface_edge);
    free(dedge_vtx_idx);
    free(dedge_vtx);
    free(dedge_face);
    free(dedge_group_idx);
    free(dedge_group);
  }

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
