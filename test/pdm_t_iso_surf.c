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

  // return PDM_ABS(x) + PDM_ABS(y) + PDM_ABS(z) - 0.5;

  // double v1 = x * x + y * y + z * z - 0.125;
  // double v2 = (x-0.2) * (x-0.2) + (y-0.2) * (y-0.2) + (z-0.3) * (z-0.3) - 0.02;
  // return PDM_MIN(v1, v2);

  // return cos(6*x * y * z) - 0.999;
}
// coordsX * coordsX + coordsY * coordsY + coordsZ * coordsZ - 0.125
// (coordsX-0.2) * (coordsX-0.2) + (coordsY-0.2) * (coordsY-0.2) + (coordsZ-0.3) * (coordsZ-0.3) - 0.07

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

  // *df_dx = PDM_SIGN(x);
  // *df_dy = PDM_SIGN(y);
  // *df_dz = PDM_SIGN(z);

  // double v1 = x * x + y * y + z * z - 0.125;
  // double v2 = (x-0.2) * (x-0.2) + (y-0.2) * (y-0.2) + (z-0.3) * (z-0.3) - 0.02;

  // if (v1 < v2) {
  //   *df_dx = 2*x;
  //   *df_dy = 2*y;
  //   *df_dz = 2*z;
  // } else {
  //   *df_dx = 2*(x-0.2);
  //   *df_dy = 2*(y-0.2);
  //   *df_dz = 2*(z-0.3);
  // };

  // double s = -6*sin(6*x*y*z);
  // *df_dx = y*z*s;
  // *df_dy = x*z*s;
  // *df_dz = x*y*s;
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
             &length,
             &elt_type);

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
   *  Create distributed cube
   */

  PDM_g_num_t *distrib_cell = NULL;
  PDM_g_num_t *distrib_face = NULL;
  PDM_g_num_t *distrib_edge = NULL;
  PDM_g_num_t *vtx_distrib  = NULL;

  PDM_dcube_nodal_t          *dcube   = NULL;
  PDM_dmesh_nodal_to_dmesh_t *dmntodm = NULL;

  PDM_g_num_t  ng_cell = 0;
  PDM_g_num_t  ng_face = 0;
  PDM_g_num_t  ng_vtx  = 0;
  int          dn_cell = 0;
  int          dn_face = 0;
  int          dn_edge = 0;
  int          dn_vtx  = 0;
  double      *dvtx_coord     = NULL;
  int         *dcell_face_idx = NULL;
  PDM_g_num_t *dcell_face     = NULL;
  int         *dface_edge_idx = NULL;
  PDM_g_num_t *dface_edge     = NULL;
  int         *dedge_vtx_idx  = NULL;
  PDM_g_num_t *dedge_vtx      = NULL;
  int          n_face_group = 0;
  PDM_g_num_t *dface_cell      = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  if (elt_type < PDM_MESH_NODAL_POLY_3D) {
    dcube = PDM_dcube_nodal_gen_create (comm,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        length,
                                        -0.5*length,//-0.1,
                                        -0.5*length,//-0.2,
                                        -0.5*length,//0.07,
                                        elt_type,
                                        1,
                                        PDM_OWNERSHIP_KEEP);

    PDM_dcube_nodal_gen_build (dcube);

    PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
    PDM_dmesh_nodal_generate_distribution(dmn);

    vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

    if (1) {
      double noise = 0.2 * length / (double) (n_vtx_seg - 1);
      for (int i = 0; i < dn_vtx; i++) {
        for (int j = 0; j < 3; j++) {
          dvtx_coord[3*i+j] += noise * (rand() / (double) RAND_MAX - 0.5);
        }
        double y = dvtx_coord[3*i+1];
        dvtx_coord[3*i] += 0.15*length*sin(7*y);
        double x = dvtx_coord[3*i];
        dvtx_coord[3*i+2] += 0.2*length*cos(8*x);
      }
    }

    if(1 == 1) {
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
    }

    dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

    PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

    PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

    PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

    PDM_dmesh_t* dmesh = NULL;
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
  }

  else {
    // polyvol_gen
    PDM_poly_vol_gen (comm,
                      -0.1,//-0.5,
                      0.1,//-0.5,
                      0.13,//-0.5,
                      length,
                      length,
                      length,
                      n_vtx_seg,
                      n_vtx_seg,
                      n_vtx_seg,
                      1,
                      0,
                      &ng_cell,
                      &ng_face,
                      &ng_vtx,
                      &n_face_group,
                      &dn_cell,
                      &dn_face,
                      &dn_vtx,
                      &dcell_face_idx,
                      &dcell_face,
                      &dface_cell,
                      &dface_edge_idx,
                      &dface_vtx,
                      &dvtx_coord,
                      &dface_group_idx,
                      &dface_group);

    printf("TO DO: build edges\n");
    abort();
  }


  // Compute dfield and gradient field
  double *dfield          = (double *) malloc(     dn_vtx * sizeof(double));
  double *dgradient_field = (double *) malloc( 3 * dn_vtx * sizeof(double));

  for(int i = 0; i < dn_vtx; ++i) {

    double x1 = dvtx_coord[3*i  ];
    double y1 = dvtx_coord[3*i+1];
    double z1 = dvtx_coord[3*i+2];
    dfield[i] = _unit_sphere(x1, y1, z1);

    _unit_sphere_gradient(x1, y1, z1,
                          &dgradient_field[3*i],
                          &dgradient_field[3*i+1],
                          &dgradient_field[3*i+2]);
  }


  PDM_iso_surface_t* isos = PDM_iso_surface_create(3, PDM_ISO_SURFACE_KIND_FIELD, 1, PDM_OWNERSHIP_KEEP, comm);
  // PDM_iso_surface_t* isos = PDM_iso_surface_create(3, PDM_ISO_SURFACE_KIND_PLANE, 1, PDM_OWNERSHIP_KEEP, comm);

  PDM_iso_surface_plane_equation_set(isos, 1., 0., 0., -0);
  // PDM_iso_surface_plane_equation_set(isos, 1., 0.5, 0.25, -0.0234);

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
  PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_VERTEX, vtx_distrib);

  PDM_iso_surface_dvtx_coord_set (isos, dvtx_coord     );
  PDM_iso_surface_dfield_set     (isos, dfield         );
  PDM_iso_surface_dgrad_field_set(isos, dgradient_field);

  PDM_iso_surface_compute(isos);

  PDM_iso_surface_free(isos);

  free(dfield);
  free(dgradient_field);


  if (elt_type < PDM_MESH_NODAL_POLY_3D) {
    PDM_dmesh_nodal_to_dmesh_free(dmntodm);
    PDM_dcube_nodal_gen_free(dcube);
  } else {
    free(distrib_cell);
    free(distrib_face);
    free(distrib_edge);
    free(vtx_distrib);
    free(dvtx_coord);
    free(dcell_face_idx);
    free(dcell_face);
    free(dface_edge_idx);
    free(dface_edge);
    free(dedge_vtx_idx);
    free(dedge_vtx);
    free(dface_cell);
    free(dface_vtx);
    free(dface_group_idx);
    free(dface_group);
  }

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

}
