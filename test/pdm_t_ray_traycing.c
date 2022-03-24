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
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"
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
#include "pdm_plane.h"
#include "pdm_triangle.h" 


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
void
_generate_surface_mesh
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         nu,
 const PDM_g_num_t         nv,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
 const PDM_split_dual_t    part_method,
 const int                 n_part,
       PDM_dmesh_nodal_t **_dmn,
       PDM_multipart_t   **_mpart
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  assert(nu > 1);
  assert(nv > 1);

  /*
   *  Vertices
   */
  PDM_g_num_t gn_vtx = (nu - 1) * nv + 2;

  PDM_g_num_t *distrib_vtx = PDM_compute_uniform_entity_distribution(comm,
                                                                     gn_vtx);

  const double step_u = 2*PDM_PI / (double) (nu - 2);
  const double step_v =   PDM_PI / (double) (nv + 1);

  int dn_vtx = (int) (distrib_vtx[i_rank+1] - distrib_vtx[i_rank]);
  double *dvtx_coord = (double *) malloc(sizeof(double) * dn_vtx * 3);

  for (int i = 0; i < dn_vtx; i++) {

    PDM_g_num_t g = distrib_vtx[i_rank] + i;

    if (g == gn_vtx-1) {
      // north pole
      dvtx_coord[3*i  ] = x_center;
      dvtx_coord[3*i+1] = y_center;
      dvtx_coord[3*i+2] = z_center + radius;

    } else if (g == gn_vtx-2) {
      // south pole
      dvtx_coord[3*i  ] = x_center;
      dvtx_coord[3*i+1] = y_center;
      dvtx_coord[3*i+2] = z_center - radius;

    } else {

      PDM_g_num_t jj = g / (nu - 1);
      PDM_g_num_t ii = g % (nu - 1);

      double v = -0.5*PDM_PI + (jj+1)*step_v;
      double u = ii*step_u;
      double c = cos(v);

      dvtx_coord[3*i  ] = x_center + radius * c * cos(u);
      dvtx_coord[3*i+1] = y_center + radius * c * sin(u);
      dvtx_coord[3*i+2] = z_center + radius * sin(v);

    }

  }

  free(distrib_vtx);



  /*
   *  Faces
   */
  PDM_g_num_t gn_face = 2*((nu - 1)*(nv - 1) + nu - 1);

  PDM_g_num_t *distrib_face = PDM_compute_uniform_entity_distribution(comm,
                                                                      gn_face);

  int          dn_face       = (int) (distrib_face[i_rank+1] - distrib_face[i_rank]);
  int         *dface_vtx_idx = (int *)         malloc(sizeof(int)         * (dn_face + 1));
  PDM_g_num_t *dface_vtx     = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dn_face * 3);

  dface_vtx_idx[0] = 0;
  for (int i = 0; i < dn_face; i++) {

    PDM_g_num_t g = distrib_face[i_rank] + i;

    dface_vtx_idx[i+1] = dface_vtx_idx[i] + 3;
    PDM_g_num_t *_fv = dface_vtx + dface_vtx_idx[i];


    if (g >= 2*(nu - 1)*(nv - 1) + nu - 1) {

      // north pole cap
      PDM_g_num_t ii = g - (2*(nu - 1)*(nv - 1) + nu - 1);

      _fv[0] = 1 + ii            + (nu-1)*(nv-1);
      _fv[1] = 1 + (ii+1)%(nu-1) + (nu-1)*(nv-1);
      _fv[2] = gn_vtx;
    }

    else if (g >= 2*(nu - 1)*(nv - 1)) {

      // south pole cap
      PDM_g_num_t ii = g - 2*(nu - 1)*(nv - 1);

      _fv[0] = 1 + ii;
      _fv[1] = gn_vtx - 1;
      _fv[2] = 1 + (ii+1)%(nu-1);
    }

    else {

      if (g >= (nu - 1)*(nv - 1)) {
        g -= (nu - 1)*(nv - 1);

        PDM_g_num_t jj = g / (nu - 1);
        PDM_g_num_t ii = g % (nu - 1);

        _fv[0] = 1 + (ii+1)%(nu-1) + (nu-1)*jj;
        _fv[1] = 1 + (ii+1)%(nu-1) + (nu-1)*(jj+1);
        _fv[2] = 1 + ii            + (nu-1)*(jj+1);
      }

      else {

        PDM_g_num_t jj = g / (nu - 1);
        PDM_g_num_t ii = g % (nu - 1);

        _fv[0] = 1 + ii +            (nu-1)*jj;
        _fv[1] = 1 + (ii+1)%(nu-1) + (nu-1)*jj;
        _fv[2] = 1 + ii            + (nu-1)*(jj+1);

      }


    }

  }

  free(distrib_face);


  /*
   *  Create dmesh nodal
   */
  PDM_dmesh_nodal_t *dmn = PDM_DMesh_nodal_create(comm,
                                                  2,
                                                  gn_vtx,
                                                  0,
                                                  gn_face,
                                                  0);

  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  dmn->surfacic->n_g_elmts = gn_face;
  int id_section = PDM_DMesh_nodal_elmts_section_add(dmn->surfacic,
                                                     PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_elmts_section_std_set(dmn->surfacic,
                                        id_section,
                                        dn_face,
                                        dface_vtx,
                                        PDM_OWNERSHIP_KEEP);

  int n_group = 0;
  int *dgroup_elt_idx = (int *) malloc(sizeof(int) * (n_group + 1));
  dgroup_elt_idx[0] = 0;
  PDM_g_num_t *dgroup_elt = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dgroup_elt_idx[n_group]);
  PDM_DMesh_nodal_elmts_group_set(dmn->surfacic,
                                  n_group,
                                  dgroup_elt_idx,
                                  dgroup_elt,
                                  PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_generate_distribution(dmn);

  /*
   * Partitionnement
   */
  int n_zone = 1;
  // int n_part_zones = {n_part};
  int *n_part_zones = (int *) malloc(sizeof(int) * n_zone);
  n_part_zones[0] = n_part;

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

  *_mpart = mpart;
  *_dmn   = dmn;
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
                                                         -0.5*length,
                                                         -0.5*length,
                                                         -0.5*length,
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
   *  Generate surface mesh
   */
  if (i_rank == 0) {
    printf("-- Generate surface mesh\n");
    fflush(stdout);
  }
  PDM_g_num_t nu = 2*n_vtx_seg;
  PDM_g_num_t nv = n_vtx_seg;

  PDM_dmesh_nodal_t     *dmn_surf   = NULL;
  PDM_multipart_t       *mpart_surf = NULL;
  _generate_surface_mesh (comm,
                          nu,
                          nv,
                          0.,
                          0.,
                          0.,
                          0.3*length,
                          part_method,
                          n_part,
                          &dmn_surf,
                          &mpart_surf);

  int          *surf_pn_vtx         = (int *)          malloc(sizeof(int)           * n_part);
  int          *surf_pn_face        = (int *)          malloc(sizeof(int)           * n_part);
  double      **surf_pvtx_coord     = (double **)      malloc(sizeof(double *)      * n_part);
  PDM_g_num_t **surf_pvtx_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **surf_pface_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    surf_pn_vtx[i_part] = PDM_multipart_part_vtx_coord_get(mpart_surf,
                                                           0,
                                                           i_part,
                                                           &surf_pvtx_coord[i_part],
                                                           PDM_OWNERSHIP_KEEP);

    PDM_multipart_part_ln_to_gn_get(mpart_surf,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_VERTEX,
                                    &surf_pvtx_ln_to_gn[i_part],
                                    PDM_OWNERSHIP_KEEP);

    surf_pn_face[i_part] = PDM_multipart_part_ln_to_gn_get(mpart_surf,
                                                           0,
                                                           i_part,
                                                           PDM_MESH_ENTITY_FACE,
                                                           &surf_pface_ln_to_gn[i_part],
                                                           PDM_OWNERSHIP_KEEP);
  }

  PDM_part_mesh_nodal_elmts_t* pmne_surf = PDM_dmesh_nodal_to_part_mesh_nodal_elmts(dmn_surf,
                                                                                    PDM_GEOMETRY_KIND_SURFACIC,
                                                                                    n_part,
                                                                                    surf_pn_vtx,
                                                                                    surf_pvtx_ln_to_gn,
                                                                                    surf_pn_face,
                                                                                    surf_pface_ln_to_gn,
                                                                                    NULL);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int id_section = 0;
    PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get(pmne_surf, id_section);
    int         *elmt_vtx                 = NULL;
    int         *parent_num               = NULL;
    PDM_g_num_t *numabs                   = NULL;
    PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
    PDM_part_mesh_nodal_elmts_block_std_get(pmne_surf, id_section, i_part, &elmt_vtx, &numabs, &parent_num, &parent_entitity_ln_to_gn);

    char filename[999];
    sprintf(filename, "out_surfacic_%i_%i.vtk", i_part, i_rank);
    PDM_vtk_write_std_elements(filename,
                               surf_pn_vtx[i_part],
                               surf_pvtx_coord[i_part],
                               surf_pvtx_ln_to_gn[i_part],
                               t_elt,
                               surf_pn_face[i_part],
                               elmt_vtx,
                               surf_pface_ln_to_gn[i_part],
                               0,
                               NULL,
                               NULL);
  }


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
    part_n_elt[i_part] = surf_pn_face[i_part];
    part_elt_extents[i_part] = (double *) malloc( 6 * part_n_elt[i_part] * sizeof(double));

    // double      *vtx_coord    = surf_pvtx_coord[i_part];
    // int         *face_vtx_idx = surf_pface_vtx_idx[i_part];
    // int         *face_vtx     = surf_pface_vtx[i_part];
    // PDM_g_num_t *face_g_num   = surf_pface_ln_to_gn[i_part];
    int id_section = 0;
    PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get(pmne_surf, id_section);
    int         *face_vtx                 = NULL;
    int         *parent_num               = NULL;
    PDM_g_num_t *face_g_num               = NULL;
    PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
    PDM_part_mesh_nodal_elmts_block_std_get(pmne_surf,
                                            id_section,
                                            i_part,
                                            &face_vtx,
                                            &face_g_num,
                                            &parent_num,
                                            &parent_entitity_ln_to_gn);
    part_elt_g_num[i_part] = face_g_num;

    int stride = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, 1);

    double *_extents = part_elt_extents[i_part];
    for (int j = 0; j < part_n_elt[i_part]; j++) {

      for (int k1 = 0; k1 < 3; k1++) {
        _extents[  k1] = DBL_MAX;
        _extents[3+k1] = -DBL_MAX;
      }

      for (int k = stride*j; k < stride*(j+1); k++) {
        int i_vtx = face_vtx[k] - 1;
        double *_coords = (double *) surf_pvtx_coord[i_part] + 3 * i_vtx;

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
        _extents[k1]   -= delta;
        _extents[3+k1] += delta;
      }
      _extents += 6;
    }


    if (1) {
      char filename[999];
      sprintf(filename, "part_elt_extents_%i_%i.vtk", i_part, i_rank);
      PDM_vtk_write_boxes(filename,
                          part_n_elt[i_part],
                          part_elt_extents[i_part],
                          face_g_num);
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

  log_trace("global_extents = %f %f %f   %f %f %f\n",
            global_extents[0], global_extents[1], global_extents[2],
            global_extents[3], global_extents[4], global_extents[5]);

  PDM_dbbtree_t *dbbt = PDM_dbbtree_create (comm, 3, global_extents);


  PDM_box_set_t  *surf_mesh_boxes = PDM_dbbtree_boxes_set (dbbt,
                                                           n_part_mesh,
                                                           part_n_elt,
                                      (const double **)    part_elt_extents,
                                  (const PDM_g_num_t **)   part_elt_g_num);


  // int          n_line = 0;
  // PDM_g_num_t *line_g_num = malloc(6 * n_line * sizeof(PDM_g_num_t));
  // double      *line_coord = malloc(6 * n_line * sizeof(double     ));

  /*
   *  Construction des lignes
   */
  int          *pn_selected_vtx      = malloc(n_part * sizeof(int        ));
  PDM_g_num_t **p_selected_vtx_g_num = malloc(n_part * sizeof(PDM_g_num_t *));
  double      **p_selected_vtx_coord = malloc(n_part * sizeof(double *));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    // int id_section = 0;
    // PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get(pmne_vol, id_section);
    // int         *elmt_vtx                 = NULL;
    // int         *parent_num               = NULL;
    // PDM_g_num_t *numabs                   = NULL;
    // PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
    // PDM_part_mesh_nodal_elmts_block_std_get(pmne_vol, id_section, i_part, &elmt_vtx, &numabs, &parent_num, &parent_entitity_ln_to_gn);

    pn_selected_vtx     [i_part] = 0;
    p_selected_vtx_g_num[i_part] = malloc(pn_vtx[i_part] * sizeof(PDM_g_num_t));
    p_selected_vtx_coord[i_part] = malloc(3 * pn_vtx[i_part] * sizeof(double));

    for(int i_vtx = 0; i_vtx < pn_vtx[i_part]; ++i_vtx) {

      int is_inside = 1;
      for (int j = 0; j < 3; j++) {
        if (pvtx_coord[i_part][3*i_vtx+j] < global_extents[j] ||
            pvtx_coord[i_part][3*i_vtx+j] > global_extents[j+3]) {
          is_inside = 0;
        }
      }

      if(is_inside == 1) {
        memcpy(p_selected_vtx_coord[i_part] + 3*pn_selected_vtx[i_part],
               pvtx_coord[i_part] + 3*i_vtx,
               sizeof(double) * 3);
        p_selected_vtx_g_num[i_part][pn_selected_vtx[i_part]] = pvtx_ln_to_gn[i_part][i_vtx];
        pn_selected_vtx[i_part]++;
      }

    }

    p_selected_vtx_g_num[i_part] = realloc(p_selected_vtx_g_num[i_part], pn_selected_vtx[i_part] * sizeof(PDM_g_num_t));

  }

  if(1 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_long(p_selected_vtx_g_num[i_part], pn_selected_vtx[i_part], "p_selected_vtx_g_num ::") ;
    }
  }


  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                                      p_selected_vtx_g_num,
                                                      NULL,
                                                      pn_selected_vtx,
                                                      n_part,
                                                      comm);

  PDM_g_num_t *distrib_selected_vtx_idx = PDM_part_to_block_distrib_index_get(ptb);
  PDM_UNUSED(distrib_selected_vtx_idx);
  int dn_vtx_selected = PDM_part_to_block_n_elt_block_get (ptb);
  PDM_g_num_t *block_g_num = PDM_part_to_block_block_gnum_get (ptb);

  PDM_log_trace_array_long(block_g_num, dn_vtx_selected, "block_g_num : ");

  double *blk_vtx_coord = NULL;
  PDM_part_to_block_exch (ptb,
                          3 * sizeof(double),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) p_selected_vtx_coord,
                          NULL,
                (void **) &blk_vtx_coord);
  for (int i_part = 0; i_part < n_part; i_part++) {
    free(p_selected_vtx_coord[i_part]);
  }
  free(p_selected_vtx_coord);


  if (1) {
    char filename[999];
    sprintf(filename, "blk_vtx_coord_%2.2d.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              dn_vtx_selected,
                              blk_vtx_coord,
                              NULL,
                              NULL);
  }


  /*
   *  On a selectionné pas mal de vertex maintenant on :
   *     -calcul les rayons
   */
  int    n_rayons = 5;
  double rayon_bb = 1.2*sqrt(pow(global_extents[3]-global_extents[0],2)
                            +pow(global_extents[4]-global_extents[1],2)
                            +pow(global_extents[5]-global_extents[2],2));
  int n_line=n_rayons*dn_vtx_selected;
  double  *line_coord = malloc(6 * n_line * sizeof(double));
  PDM_g_num_t *line_g_num = malloc(n_line * sizeof(PDM_g_num_t));
  PDM_g_num_t *distrib = PDM_part_to_block_distrib_index_get(ptb);
  double alpha;
  double phi   = 0.;
  double theta = 0.;
  for (int i_vtx = 0; i_vtx < dn_vtx_selected; ++i_vtx) {

    for (int j_rayon = 0; j_rayon < n_rayons; ++j_rayon) {

      line_g_num[n_rayons*i_vtx+j_rayon] = distrib[i_rank]+n_rayons*i_vtx+j_rayon+1;
      
      if (j_rayon==1){
        alpha = (double) rand() / (double) RAND_MAX;
        phi=alpha*PDM_PI;
        alpha = (double) rand() / (double) RAND_MAX;
        theta=alpha*2.*PDM_PI;
      } 
      else if (j_rayon==2){
        phi=phi+PDM_PI*0.5;
      } 
      else if (j_rayon==3){
        theta=theta-PDM_PI*0.5;
        phi=phi-PDM_PI*0.5;
      } 
      else if (j_rayon==4){
        phi=phi+PDM_PI*0.5;
      } 
      else if (j_rayon==5){
        theta=theta+3.*PDM_PI*0.5;
        phi=phi-PDM_PI*0.5;
      } 
      else if (j_rayon==6){
        theta=theta-PDM_PI*0.5;
      } 
      else if ((j_rayon>=7) && (j_rayon<=14)){
        theta=(j_rayon-7)*PDM_PI*0.25;
        phi=PDM_PI*0.25;
      } 
      else if ((j_rayon>=15) && (j_rayon<=18)){
        theta=(2*(j_rayon-15)+1)*PDM_PI*0.25;
        phi=PDM_PI*0.5;
      } 
      else if ((j_rayon>=19) && (j_rayon<=26)){
        theta=(j_rayon-19)*PDM_PI*0.25;
        phi=3.*PDM_PI*0.5;
      } 
      else
      { 
        alpha = (double) rand() / (double) RAND_MAX;
        phi=alpha*PDM_PI;
        alpha = (double) rand() / (double) RAND_MAX;
        theta=alpha*2.*PDM_PI;
      } 
      
      line_coord[0+6*j_rayon+6*n_rayons*i_vtx]=blk_vtx_coord[0+3*i_vtx]; 
      line_coord[1+6*j_rayon+6*n_rayons*i_vtx]=blk_vtx_coord[1+3*i_vtx];
      line_coord[2+6*j_rayon+6*n_rayons*i_vtx]=blk_vtx_coord[2+3*i_vtx];
      //  Calcul des coordonnees du point de controle
      line_coord[3+6*j_rayon+6*n_rayons*i_vtx]=line_coord[0+6*j_rayon+6*n_rayons*i_vtx]+rayon_bb*cos(theta)*sin(phi);
      line_coord[4+6*j_rayon+6*n_rayons*i_vtx]=line_coord[1+6*j_rayon+6*n_rayons*i_vtx]+rayon_bb*sin(theta)*sin(phi);
      line_coord[5+6*j_rayon+6*n_rayons*i_vtx]=line_coord[2+6*j_rayon+6*n_rayons*i_vtx]+rayon_bb*cos(phi);

    } 
  }


  int         *box_idx   = NULL;
  PDM_g_num_t *box_g_num = NULL;
  PDM_dbbtree_lines_intersect_boxes(dbbt,
                                    n_line,
                                    line_g_num,
                                    line_coord,
                                    &box_idx,
                                    &box_g_num);

  double *weight = malloc(sizeof(double)*box_idx[n_line]);

  for (int i = 0; i < box_idx[n_line]; ++i) {
    weight[i]=1.;
  } 

  /*
   * @Bastien :
   *   L'idée ici est d'avoir pour chaque boite le nombre de rayons qui la traverse
   *   Si nb_rayons = 0, on degage la boite
   */
  // int          n_select_boxes   = 0;
  // double      *box_center_coord = NULL;
  // // PDM_g_num_t *box_g_num        = NULL;
  // int         *box_ray_n        = NULL;
  // PDM_part_to_block_t* ptb_equi = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
  //                                                               PDM_PART_TO_BLOCK_POST_CLEANUP,
  //                                                               1.,
  //                                                               PDM_PART_GEOM_HILBERT,
  //                                                               &box_center_coord,
  //                                                               &box_g_num,
  //                                                               &box_ray_n, // =  Weight
  //                                                               &n_select_boxes,
  //                                                               1,
  //                                                               comm);

  // int          n_parent         = PDM_part_to_block_n_elt_block_get     (ptb_equi);
  // PDM_g_num_t *parent_gnum      = PDM_part_to_block_block_gnum_get      (ptb_equi);
  // PDM_g_num_t *distrib_equi_box = PDM_part_to_block_distrib_index_get   (ptb_equi);
  // Il faut s'assurer qu'on envoie que des boites qui ont au moins un rayons
  // int         *equi_box_count   = PDM_part_to_block_block_gnum_count_get(ptb_equi);

  /*
   *  Pour la localisation : on peut échanger le type d'élement puis on fait une permutation
   *  sur le parent_gnum (avec un idx par type )
   */


  /*
   * Echange de l'information box->rayon dans le nouveau partitionnement
   *   See : paradigm/test/pdm_t_extract_part_from_part.c ( A adapter car on a deja lesconectivités en gnum )
   *   PDM_pconnectivty_to_pconnectivty(pbox_ray, pbox_ray_idx) --> ray_ln_to_gn + box_ray
   *   part1 = box
   *   part2 = ray
   *   part1_to_part2 = box_ray_idx + box_ray_gnum ( et donc écahnge reverse )
   *   Puis on échange aussi les rayons (c'est le même part_to_part)
   *   Il faudrait découper PDM_pconnectivty_to_pconnectivty pour pouvoir garder le part_to_part
   *      - Creation du protocol
   *      - Echange spécifique a des connectivités (pour les passés en local )
   *      - Finalisation
   *
   *
   *  DOnc a ce stade, on a une connectivité local box_ray et les info de rayons locals (pas de duplication !!!!)
   */


  /*
   * On peut faire je pense avec le ptb_equi mais :
   *      - pour la localisation on peut pas faire le tri indirect pour choper des sections
   *      - Plus dure de renvoyer les résultats après
   * Nouvaux part_to_part entre les box redistribué et le partitionnement initiale
   * Donc part1 = implicite ordering of current block
   *      part2 = user ordering
   *      part1_to_part2 = parent_gnum (si mesh location il est permuté pour avoir des block d'elements en entrée )
   * Echange des connectivités part part_to_part reverse
   */

  /*
   *  On a en local tout pour faire la localisation / ray_traycing
   */

  /*
   * Renvoie des résultats par part_to_part
   *   --> Direct
   */



  PDM_part_to_block_t *ptb_face = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                           1.,
                                                          &box_g_num,
                                                          &weight,
                                                          &box_idx[n_line],
                                                           1,
                                                           comm);

  PDM_g_num_t *distrib_faces = PDM_part_to_block_distrib_index_get(ptb_face);

  PDM_g_num_t *part_data_g_num = malloc(sizeof(PDM_g_num_t)*box_idx[n_line]);
  int *part_stride = malloc(sizeof(int)*box_idx[n_line]);
  double* part_data_coord = malloc(6*sizeof(double)*box_idx[n_line]);


  for (int i = 0; i < n_line; ++i) {
    for (int j = box_idx[i] ; j < box_idx[i+1] ; ++j) {
      part_data_g_num[j]=line_g_num[i];
      part_stride[j]=1;
      for (int k = 0; k < 6; ++k) {
        part_data_coord[6*j+k]=line_coord[6*i+k];  
      } 
    } 
  } 

  PDM_g_num_t *dface_line = NULL;
  int *dface_line_n =NULL;

  PDM_part_to_block_exch (ptb_face,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                         &part_stride,
              (void **)  &part_data_g_num,
                         &dface_line_n,
              (void **)  &dface_line);
  free(dface_line_n);
  double *dface_line_coord =NULL;

  PDM_part_to_block_exch (ptb_face,
                          6 * sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                         &part_stride,
              (void **)  &part_data_coord,
                         &dface_line_n,
              (void **)  &dface_line_coord);

  int dn_face=PDM_part_to_block_n_elt_block_get(ptb_face);
  int* dface_line_idx = PDM_array_new_idx_from_sizes_int(dface_line_n,dn_face);

  PDM_log_trace_connectivity_long(dface_line_idx,dface_line,dn_face,"dface_line :");

  free(dface_line_n);

  /*******************************************************************
    *  Transfer element coords from parts to blocks
    *******************************************************************/

  PDM_g_num_t l_max_elt_g_num = 0;

  for (int i_part = 0; i_part < n_part; i_part++) {

    int          id_section = 0;
    int         *elmt_vtx                 = NULL;
    int         *parent_num               = NULL;
    PDM_g_num_t *numabs                   = NULL;
    PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
    PDM_part_mesh_nodal_elmts_block_std_get(pmne_surf, id_section, i_part, &elmt_vtx, &numabs, &parent_num, &parent_entitity_ln_to_gn);
    for (int j_part = 0 ; j_part<surf_pn_face[i_part] ; j_part++){
      l_max_elt_g_num=PDM_MAX(l_max_elt_g_num,numabs[j_part]);
    } 
  } 

  PDM_g_num_t g_max_elt_g_num = 0;


  PDM_MPI_Allreduce (&l_max_elt_g_num,  &g_max_elt_g_num,   1,
                    PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);

  PDM_g_num_t *_block_elt_distrib_idx = distrib_faces;


  if (distrib_faces[n_rank] < g_max_elt_g_num) {
    _block_elt_distrib_idx = malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
    for (int i = 0; i < n_rank; i++) {
      _block_elt_distrib_idx[i] = distrib_faces[i];
    }
    _block_elt_distrib_idx[n_rank] = g_max_elt_g_num;
  }

  PDM_part_to_block_t *ptb_elt = PDM_part_to_block_create_from_distrib (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                        PDM_PART_TO_BLOCK_POST_MERGE,
                                                                        1.,
                                                                        (PDM_g_num_t **) surf_pface_ln_to_gn,
                                                                        _block_elt_distrib_idx,
                                                                        surf_pn_face,
                                                                        n_part_mesh,
                                                                        comm);

  int **part_elt_vtx_n = malloc (sizeof(int *) * n_part_mesh);
  double **part_elt_vtx_coord = malloc (sizeof(double *) * n_part_mesh);
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {

    int id_section = 0;
    PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get(pmne_surf, id_section);
    int         *part_face_vtx                 = NULL;
    int         *parent_num               = NULL;
    PDM_g_num_t *face_g_num               = NULL;
    PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
    PDM_part_mesh_nodal_elmts_block_std_get(pmne_surf,
                                            id_section,
                                            ipart,
                                            &part_face_vtx,
                                            &face_g_num,
                                            &parent_num,
                                            &parent_entitity_ln_to_gn);

    int stride = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, 1);

    const double *part_vtx = surf_pvtx_coord[ipart];

    part_elt_vtx_n[ipart] = malloc (sizeof(int) * surf_pn_face[ipart]);
    part_elt_vtx_coord[ipart] = malloc (sizeof(double) * stride * surf_pn_face[ipart] * 3);

    for (int i = 0; i < surf_pn_face[ipart]; i++) {
      part_elt_vtx_n[ipart][i] = stride;//3 * face_vtx_n;
      for (int j = stride*i; j < stride*(i+1); j++) {
        int ivtx = part_face_vtx[j] - 1;
        for (int k = 0; k < 3; k++) {
          part_elt_vtx_coord[ipart][3*j + k] = part_vtx[3*ivtx + k];
        }
      }
    }
  }

  int *block_elt_vtx_n = NULL;
  double *block_elt_vtx_coord = NULL;
  PDM_part_to_block_exch (ptb_elt,
                          3*sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          -1,
                          part_elt_vtx_n,
                          (void **) part_elt_vtx_coord,
                          &block_elt_vtx_n,
                          (void **) &block_elt_vtx_coord);
  for (int ipart = 0; ipart < n_part; ipart++) {
    free (part_elt_vtx_n[ipart]);
    free (part_elt_vtx_coord[ipart]);
  }
  free (part_elt_vtx_n);
  free (part_elt_vtx_coord);


  int* block_elt_vtx_idx = PDM_array_new_idx_from_sizes_int(block_elt_vtx_n,dn_face);


  //dn_face : nb faces dans un block 
  //graph face rayon dface_line (donne les n absolus des rayons candidats)
  //dface_line_idx index de parcour de dface_line
  // définition des faces donné par block_elt_vtx_n
  // block_elt_vtx_coord coordonnés des vertex des faces 

  double epsilon = 1e-12;

  int  *tag_face = malloc(dface_line_idx[dn_face] * sizeof(int));

  for (int i_face = 0; i_face < dn_face; ++i_face) {
    double face_normal[3];
    double* face_coord=block_elt_vtx_coord+3*block_elt_vtx_idx[i_face];

    PDM_plane_normal(block_elt_vtx_n[i_face],face_coord,face_normal);
    for (int j_rayon = dface_line_idx[i_face]; j_rayon < dface_line_idx[i_face+1]; ++j_rayon) {

      double* rayon_coord=dface_line_coord+6*j_rayon;
      int stat = -1;
      double intersection[3];
      double vec[3] = {
        face_coord[0] - rayon_coord[0],
        face_coord[1] - rayon_coord[1],
        face_coord[2] - rayon_coord[2]
      };

      double rayon_vect[3] = {
        rayon_coord[3] - rayon_coord[0],
        rayon_coord[4] - rayon_coord[1],
        rayon_coord[5] - rayon_coord[2]
      };

      double denom = PDM_DOT_PRODUCT(rayon_vect, face_normal);
      double numer = PDM_DOT_PRODUCT(vec, face_normal);

      if (PDM_ABS(denom) < epsilon) {   // Epsilon is here to avoid division by 0
        if (PDM_ABS(numer) < epsilon) { // Le edge contenu dans le plan de la face
          stat = 2;
          memcpy(intersection, rayon_coord, sizeof(double) * 3);
        } else {// Edge parallèle mais pas dans le même plan
          stat = 0;
        }
      } else {  // Ca intersecte nickel
        double t = numer/denom;

        if ((t >=0)&&(t<1)) {
          stat = 1;
          for (int j = 0; j < 3; j++) {
            intersection[j] = rayon_coord[j] + t*rayon_vect[j];
          }
        }
      }
      if(stat>0){
        double closestPoint[3]; 
        double min_dist2;
        PDM_triangle_status_t statriangle= PDM_triangle_evaluate_position(intersection,
                                                                          face_coord,
                                                                          closestPoint,
                                                                         &min_dist2,
                                                                          NULL);
        if(statriangle!=PDM_TRIANGLE_INSIDE){
          tag_face[j_rayon]=0;
        } 
        else{
          tag_face[j_rayon]=1; 
        } 
      } 
    } 
  }



  // distrib_faces[0]=-1;

  // int         *dface_rayon_idx;
  // PDM_g_num_t *dface_rayon;
  // PDM_dconnectivity_transpose(comm,
  //                             distrib,
  //                             distrib_faces,
  //                             box_idx,
  //                             box_g_num,
  //                             1,
  //                             &dface_rayon_idx,
  //                             &dface_rayon);

  
  // int dn_faces=distrib_faces[i_rank+1]-distrib_faces[i_rank]; 

  // PDM_log_trace_connectivity_long(dface_rayon_idx,dface_rayon,dn_faces,"dface_rayon :");



  // double  *tagibc = malloc(dn_vtx_selected * sizeof(double));

  // for (int i_vtx = 0; i_vtx < dn_vtx_selected; ++i_vtx) {
  //   tagibc[i_vtx]=0.;
  //   for (int j_rayon = 0; j_rayon < n_rayons; ++j_rayon) {
      


  //   } 
  // }







  free(blk_vtx_coord);
  PDM_part_to_block_free(ptb);

  if (1) {
    for (int i = 0; i < n_line; i++) {
      log_trace("line "PDM_FMT_G_NUM": ", line_g_num[i]);
      for (int j = box_idx[i]; j < box_idx[i+1]; j++) {
        log_trace(PDM_FMT_G_NUM" ", box_g_num[j]);
      }
      log_trace("\n");
    }
  }



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
  PDM_part_mesh_nodal_elmts_free(pmne_surf);

  PDM_dcube_nodal_gen_free(dcube);
  PDM_multipart_free(mpart);
  PDM_multipart_free(mpart_surf);
  PDM_DMesh_nodal_free(dmn_surf);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

}
