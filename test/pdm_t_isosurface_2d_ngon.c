#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"

#include "pdm_logging.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_vtk.h"

#include "pdm_array.h"
#include "pdm_reader_gamma.h"
#include "pdm_dcube_nodal_gen.h"

#include "pdm_isosurface.h"
#include "pdm_multipart.h"
#include "pdm_multipart_priv.h"
#include "pdm_partitioning_algorithm.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static void
_gen_mesh
(
 PDM_MPI_Comm          comm,
 const char           *filename,
 int                   n_part,
 PDM_g_num_t           n_vtx_seg,
 int                   randomize,
 PDM_Mesh_nodal_elt_t  elt_type,
 int                   generate_edges,
 PDM_multipart_t     **mpart,
 PDM_part_mesh_t      *pmesh,
 PDM_dmesh_t         **out_dmesh
 )
{
  if (PDM_Mesh_nodal_elt_dim_get(elt_type) != 2) {
    PDM_error(__FILE__, __LINE__, 0, "Element type %d is not 2D\n", elt_type);
  }

  if (n_part > 0) {
    int n_domain = 1;
    *mpart = PDM_multipart_create(n_domain,
                                  &n_part,
                                  PDM_FALSE,
                                  PDM_SPLIT_DUAL_WITH_HILBERT,
                                  PDM_PART_SIZE_HOMOGENEOUS,
                                  NULL,
                                  comm,
                                  PDM_OWNERSHIP_KEEP);
  }

  if (filename != NULL) {
    PDM_dmesh_nodal_t *dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                                          filename,
                                                          0,
                                                          0);

    if (0) {
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "isosurface_2d_ngon_surface");
    }

    if (n_part > 0) {
      PDM_multipart_dmesh_nodal_set(*mpart, 0, dmn);
      PDM_multipart_compute(*mpart);
    }
    else {
      PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1,
                                                                              comm,
                                                                              PDM_OWNERSHIP_USER);
      PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm,
                                               0,
                                               dmn);

      PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                       PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE, // not sure...
                                       PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE); // not sure...

      PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm,
                                         0,
                                         out_dmesh);

      PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);
      double *vtx_coord = PDM_DMesh_nodal_coord_get(dmn, PDM_OWNERSHIP_USER); // ugly but saving lives
      PDM_dmesh_vtx_coord_set(*out_dmesh, vtx_coord, PDM_OWNERSHIP_KEEP);
    }

    PDM_DMesh_nodal_free(dmn);
  }

  else {
    PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                          n_vtx_seg,
                                                          n_vtx_seg,
                                                          1,
                                                          1,
                                                          0,
                                                          0,
                                                          0,
                                                          elt_type,
                                                          1,
                                                          PDM_OWNERSHIP_USER);

    PDM_dcube_nodal_gen_random_factor_set(dcube, (double) randomize);

    PDM_dcube_nodal_gen_build(dcube);

    PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

    PDM_dmesh_nodal_generate_distribution(dmn);

    if (n_part > 0) {
      PDM_multipart_dmesh_nodal_set(*mpart, 0, dmn);
      PDM_multipart_compute(*mpart);
    }
    else {
      PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1,
                                                                              comm,
                                                                              PDM_OWNERSHIP_USER);
      PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm,
                                               0,
                                               dmn);

      PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                       PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                       PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);

      PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm,
                                         0,
                                         out_dmesh);

      PDM_dmesh_compute_distributions(*out_dmesh);

      PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);
      double *vtx_coord = PDM_DMesh_nodal_coord_get(dmn, PDM_OWNERSHIP_USER); // ugly but saving lives
      PDM_dmesh_vtx_coord_set(*out_dmesh, vtx_coord, PDM_OWNERSHIP_KEEP);
    }
    PDM_DMesh_nodal_free(dmn);
    PDM_dcube_nodal_gen_free(dcube);
  }



  if (pmesh != NULL) {
    assert(n_part > 0);
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i_entity = PDM_MESH_ENTITY_CELL; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
        PDM_g_num_t *entity_ln_to_gn = NULL;
        int n_entity = PDM_multipart_part_ln_to_gn_get(*mpart,
                                                       0,
                                                       i_part,
                                                       i_entity,
                                                       &entity_ln_to_gn,
                                                       PDM_OWNERSHIP_KEEP);

        PDM_part_mesh_n_entity_set(pmesh, i_part, i_entity, n_entity);
        PDM_part_mesh_entity_ln_to_gn_set(pmesh,
                                          i_part,
                                          i_entity,
                                          entity_ln_to_gn,
                                          PDM_OWNERSHIP_USER);
      }

      int *cell_face_idx = NULL;
      int *cell_face     = NULL;
      PDM_multipart_part_connectivity_get(*mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                          &cell_face_idx,
                                          &cell_face,
                                          PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_connectivity_set(pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                     cell_face,
                                     cell_face_idx,
                                     PDM_OWNERSHIP_USER);

      if (generate_edges) {
        int *face_edge_idx = NULL;
        int *face_edge     = NULL;
        PDM_multipart_part_connectivity_get(*mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                            &face_edge_idx,
                                            &face_edge,
                                            PDM_OWNERSHIP_KEEP);
        assert(face_edge != NULL);
        PDM_part_mesh_connectivity_set(pmesh,
                                       i_part,
                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                       face_edge,
                                       face_edge_idx,
                                       PDM_OWNERSHIP_USER);

        int *edge_vtx_idx = NULL;
        int *edge_vtx     = NULL;
        PDM_multipart_part_connectivity_get(*mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                            &edge_vtx_idx,
                                            &edge_vtx,
                                            PDM_OWNERSHIP_KEEP);
        assert(edge_vtx != NULL);
        PDM_part_mesh_connectivity_set(pmesh,
                                       i_part,
                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                       edge_vtx,
                                       edge_vtx_idx,
                                       PDM_OWNERSHIP_USER);
      }

      double *vtx_coord = NULL;
      PDM_multipart_part_vtx_coord_get(*mpart, 0, i_part, &vtx_coord, PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_vtx_coord_set(pmesh, i_part, vtx_coord, PDM_OWNERSHIP_USER);


      int          n_surface             = 0;
      int         *surface_face_idx      = NULL;
      int         *surface_face          = NULL;
      PDM_g_num_t *surface_face_ln_to_gn = NULL;
      PDM_multipart_group_get(*mpart,
                              0,
                              i_part,
                              PDM_MESH_ENTITY_FACE,
                              &n_surface,
                              &surface_face_idx,
                              &surface_face,
                              &surface_face_ln_to_gn,
                              PDM_OWNERSHIP_KEEP);

      PDM_part_mesh_bound_concat_set(pmesh,
                                     i_part,
                                     PDM_BOUND_TYPE_FACE,
                                     n_surface,
                                     surface_face_idx,
                                     surface_face,
                                     surface_face_ln_to_gn,
                                     PDM_OWNERSHIP_USER);
    }
  }

}


static void
_compute_iso_field
(
 int     n_vtx,
 double *vtx_coord,
 double *vtx_field
 )
{
  double ctr[3] = {0.1, 0.2, .3};
  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    double r2 = 0;
    for (int i = 0; i < 3; i++) {
      double d = vtx_coord[3*i_vtx+i] - ctr[i];
      r2 += d*d;
    }
    vtx_field[i_vtx] = 1.2*sqrt(r2);
  }
}


/**
 *
 * \brief  Usage
 *
 */
static void
_usage
(
  int exit_code
)
{
  printf
    ("\n"
     "  -h                             This message.\n\n"
     "  -n_part        <n>             Number of partitions (0 -> block-distributed).\n\n"
     "  -in            <filename>      Mesh file name (Gamma Mesh Format).\n\n"
     "  -sol           <filename>      Solution file name (Gamma Mesh Format).\n\n"
     "  -visu                          Enable exports for visualization.\n\n"
     "  -n_isovalues   <n>             Number of isovalues.\n\n"
     "  -isovalues     <v1 v2 ... vn>  Isovalues.\n\n"
     "  -elt_type      <t>             Volume element type (only for automatically generated mesh).\n\n"
     "  -randomize                     Enable mesh randomization (only for automatically generated mesh).\n\n"
     "  -n             <n>             Number of vertices per side (only for automatically generated mesh).\n\n"
     "  -use_part_mesh                 Use part_mesh structure (dmesh if n_part <= 0).\n\n"
     "  -edges                         Generate edges.\n\n"
     "  -local                         Build isosurface locally.\n\n"
     );

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 */

static void
_read_args
(
  int                    argc,
  char                 **argv,
  int                   *n_part,
  char                 **mesh_name,
  char                 **sol_name,
  int                   *visu,
  int                   *n_isovalues,
  double               **isovalues,
  PDM_Mesh_nodal_elt_t  *elt_type,
  int                   *randomize,
  PDM_g_num_t           *n_vtx_seg,
  int                   *use_part_mesh,
  int                   *generate_edges,
  int                   *local
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_part = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-in") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *mesh_name = argv[i];
      }
    }

    else if (strcmp(argv[i], "-sol") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *sol_name = argv[i];
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else if (strcmp(argv[i], "-n_isovalues") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_isovalues = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-isovalues") == 0) {
      PDM_malloc(*isovalues, *n_isovalues, double);
      for (int j = 0; j < *n_isovalues; j++) {
        i++;
        if (i >= argc)
          _usage(EXIT_FAILURE);
        else {
          (*isovalues)[j] = atof(argv[i]);
        }
      }
    }

    else if (strcmp(argv[i], "-elt_type") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-rand") == 0) {
      *randomize = 1;
    }

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg = (PDM_g_num_t) atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-use_part_mesh") == 0) {
      *use_part_mesh = 1;
    }

    else if (strcmp(argv[i], "-edges") == 0) {
      *generate_edges = 1;
    }

    else if (strcmp(argv[i], "-local") == 0) {
      *local = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */

int main
(
  int   argc,
  char *argv[]
)
{
  /*
   *  Init MPI
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  /*
   *  Read args
   */
  char                 *mesh_name      = NULL;
  char                 *sol_name       = NULL;
  int                   n_part         = 1;
  int                   visu           = 0;
  int                   n_isovalues    = 1;
  double               *isovalues      = NULL;
  PDM_Mesh_nodal_elt_t  elt_type       = PDM_MESH_NODAL_QUAD4;
  PDM_g_num_t           n_vtx_seg      = 10;
  int                   randomize      = 0;
  int                   use_part_mesh  = 0;
  int                   generate_edges = 0;
  int                   local          = 0;

  _read_args(argc,
             argv,
             &n_part,
             &mesh_name,
             &sol_name,
             &visu,
             &n_isovalues,
             &isovalues,
             &elt_type,
             &randomize,
             &n_vtx_seg,
             &use_part_mesh,
             &generate_edges,
             &local);

  if (isovalues == NULL) {
    n_isovalues = 1;
    PDM_malloc(isovalues, n_isovalues, double);
    isovalues[0] = 0.;
  }

  // Generate/read mesh
  PDM_multipart_t *mpart = NULL;
  PDM_part_mesh_t *pmesh = NULL;
  PDM_dmesh_t     *dmesh = NULL;
  if (n_part > 0) {
    if (use_part_mesh) {
      pmesh = PDM_part_mesh_create(n_part, comm);
    }
  }

  _gen_mesh(comm,
            mesh_name,
            n_part,
            n_vtx_seg,
            randomize,
            elt_type,
            generate_edges,
            &mpart,
            pmesh,
            &dmesh);


  // Compute scalar field
  double **piso_field = NULL;
  double  *diso_field = NULL;
  if (n_part > 0) {
    // Partitioned
    PDM_malloc(piso_field, n_part, double *);
    for (int i_part = 0; i_part < n_part; i_part++) {

      double *vtx_coord = NULL;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx_coord,
                                                   PDM_OWNERSHIP_KEEP);

      PDM_malloc(piso_field[i_part], n_vtx, double);
      _compute_iso_field(n_vtx, vtx_coord, piso_field[i_part]);
    }
  }
  else {
    // Block-distributed
    int dn_vtx = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_VTX);

    double *dvtx_coord = NULL;
    PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_KEEP);

    PDM_malloc(diso_field, dn_vtx, double);
    _compute_iso_field(dn_vtx, dvtx_coord, diso_field);
  }


  // Create Isosurface instance
  PDM_isosurface_t *isos = PDM_isosurface_create(comm, 2);

  if (n_part > 0) {
    PDM_extract_part_kind_t extract_kind = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    if (local) {
      extract_kind = PDM_EXTRACT_PART_KIND_LOCAL;
    }
    PDM_isosurface_redistribution_set(isos,
                                      extract_kind,
                                      PDM_SPLIT_DUAL_WITH_HILBERT);
  }


  /* Set mesh */
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  int own_dface_vtx = 0;
  if (n_part > 0) {
    // Partitioned
    if (use_part_mesh) {
      PDM_isosurface_part_mesh_set(isos, pmesh);
    }
    else {
      PDM_isosurface_n_part_set(isos, n_part);

      for (int i_part = 0; i_part < n_part; i_part++) {

        // Connectivities
        if (generate_edges) {
          int *face_edge_idx = NULL;
          int *face_edge     = NULL;
          PDM_multipart_part_connectivity_get(mpart,
                                              0,
                                              i_part,
                                              PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                              &face_edge_idx,
                                              &face_edge,
                                              PDM_OWNERSHIP_KEEP);
          PDM_isosurface_connectivity_set(isos,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                          face_edge_idx,
                                          face_edge);

          int *edge_vtx_idx = NULL;
          int *edge_vtx     = NULL;
          PDM_multipart_part_connectivity_get(mpart,
                                              0,
                                              i_part,
                                              PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                              &edge_vtx_idx,
                                              &edge_vtx,
                                              PDM_OWNERSHIP_KEEP);
          PDM_isosurface_connectivity_set(isos,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          NULL,
                                          edge_vtx);
        }
        else {
          int *face_vtx_idx = NULL;
          int *face_vtx     = NULL;
          PDM_multipart_part_connectivity_get(mpart,
                                              0,
                                              i_part,
                                              PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                              &face_vtx_idx,
                                              &face_vtx,
                                              PDM_OWNERSHIP_KEEP);

          PDM_isosurface_connectivity_set(isos,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          face_vtx_idx,
                                          face_vtx);
        }



        // Coordinates
        double *vtx_coord = NULL;
        PDM_multipart_part_vtx_coord_get(mpart,
                                         0,
                                         i_part,
                                         &vtx_coord,
                                         PDM_OWNERSHIP_KEEP);
        PDM_isosurface_vtx_coord_set(isos,
                                     i_part,
                                     vtx_coord);

        // Global IDs
        PDM_g_num_t *face_ln_to_gn = NULL;
        int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_MESH_ENTITY_FACE,
                                                     &face_ln_to_gn,
                                                     PDM_OWNERSHIP_KEEP);

        PDM_g_num_t *vtx_ln_to_gn = NULL;
        int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                    0,
                                                    i_part,
                                                    PDM_MESH_ENTITY_VTX,
                                                    &vtx_ln_to_gn,
                                                    PDM_OWNERSHIP_KEEP);

        PDM_isosurface_ln_to_gn_set(isos,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    n_face,
                                    face_ln_to_gn);

        PDM_isosurface_ln_to_gn_set(isos,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    n_vtx,
                                    vtx_ln_to_gn);

        if (generate_edges) {
          PDM_g_num_t *edge_ln_to_gn = NULL;
          int n_edge = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_MESH_ENTITY_EDGE,
                                                       &edge_ln_to_gn,
                                                       PDM_OWNERSHIP_KEEP);
          PDM_isosurface_ln_to_gn_set(isos,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      n_edge,
                                      edge_ln_to_gn);
        }
      }
    }
  }
  else {
    // Block-distributed
    if (use_part_mesh) {
      PDM_isosurface_dmesh_set(isos, dmesh);
    }
    else {
      for (int i_entity = PDM_MESH_ENTITY_FACE; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
        PDM_g_num_t *distrib = NULL;
        PDM_dmesh_distrib_get     (dmesh, i_entity, &distrib);
        PDM_isosurface_distrib_set(isos,  i_entity,  distrib);
      }

      int         *dface_edge_idx = NULL;
      PDM_g_num_t *dface_edge     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                 &dface_edge,
                                 &dface_edge_idx,
                                 PDM_OWNERSHIP_KEEP);
      int         *dedge_vtx_idx = NULL;
      PDM_g_num_t *dedge_vtx     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                 &dedge_vtx,
                                 &dedge_vtx_idx,
                                 PDM_OWNERSHIP_KEEP);
      if (generate_edges) {
        PDM_isosurface_dconnectivity_set(isos,
                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                         dface_edge_idx,
                                         dface_edge);

        PDM_isosurface_dconnectivity_set(isos,
                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                         NULL,
                                         dedge_vtx);
      }
      else {
        PDM_dmesh_connectivity_get(dmesh,
                                   PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                   &dface_vtx,
                                   &dface_vtx_idx,
                                   PDM_OWNERSHIP_KEEP);

        if (dface_vtx == NULL) {
          own_dface_vtx = 1;
          PDM_g_num_t *distrib_face = NULL;
          PDM_g_num_t *distrib_edge = NULL;
          PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE, &distrib_face);
          PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_EDGE, &distrib_edge);
          PDM_dconnectivity_dface_vtx_from_face_and_edge(comm,
                                                         distrib_face,
                                                         distrib_edge,
                                                         dface_edge_idx,
                                                         dface_edge,
                                                         dedge_vtx,
                                                         &dface_vtx);
          dface_vtx_idx = dface_edge_idx;
        }

        PDM_isosurface_dconnectivity_set(isos,
                                         PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                         dface_vtx_idx,
                                         dface_vtx);
      }

      double *dvtx_coord = NULL;
      PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_KEEP);
      PDM_isosurface_dvtx_coord_set(isos, dvtx_coord);
    }
  }


  // Plane slice
  double plane_equation [4] = {1.,0.,0.,0.5};
  double plane_isovalues[3] = {-0.30,0.,0.30};
  int iso1 = PDM_isosurface_add(isos,
                                PDM_ISO_SURFACE_KIND_PLANE,
                                3,
                                plane_isovalues);

  PDM_isosurface_equation_set(isos,
                              iso1,
                              plane_equation,
                              0);

  // Scalar field isosurface
  int iso2 = PDM_isosurface_add(isos,
                                PDM_ISO_SURFACE_KIND_FIELD,
                                n_isovalues,
                                isovalues);

  if (n_part > 0) {
    // Partitioned
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_isosurface_field_set(isos,
                               iso2,
                               i_part,
                               piso_field[i_part]);
    }
  }
  else {
    // Block-distributed
    PDM_isosurface_dfield_set(isos,
                              iso2,
                              diso_field);
  }



  /*
   *  Compute isosurface
   */
  // PDM_isosurface_compute(isos, iso1);
  // PDM_isosurface_reset  (isos, iso1);
  // PDM_isosurface_compute(isos, iso1);
  // PDM_isosurface_compute(isos, iso2);
  PDM_isosurface_compute(isos, -1);
  PDM_isosurface_reset  (isos, -1);
  PDM_isosurface_compute(isos, -1);

  /* Free memory */
  PDM_isosurface_free(isos);
  PDM_free(isovalues);

  if (n_part > 0) {
    PDM_multipart_free(mpart);

    if (use_part_mesh) {
      PDM_part_mesh_free(pmesh);
    }

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_free(piso_field[i_part]);
    }
    PDM_free(piso_field);
  }
  else {
    if (own_dface_vtx) {
      PDM_free(dface_vtx);
    }
    PDM_free(diso_field);
    PDM_dmesh_free(dmesh);
  }


  if (i_rank == 0) {
    printf("End :D\n");
    fflush(stdout);
  }



  PDM_MPI_Finalize();


  return EXIT_SUCCESS;
}
