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
#include "pdm_poly_vol_gen.h"

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
  if (PDM_Mesh_nodal_elt_dim_get(elt_type) != 3) {
    PDM_error(__FILE__, __LINE__, 0, "Element type %d is not 3D\n", elt_type);
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
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC,  "isosurface_3d_ngon_volume");
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "isosurface_3d_ngon_surface");
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
                                       PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                       PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

      PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm,
                                         0,
                                         out_dmesh);

      PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);
      double *vtx_coord = PDM_DMesh_nodal_coord_get(dmn, PDM_OWNERSHIP_USER); // ugly but saving lives
      PDM_dmesh_vtx_coord_set(*out_dmesh, vtx_coord, PDM_OWNERSHIP_KEEP);
    }

    PDM_DMesh_nodal_free(dmn);
  }

  else if (elt_type == PDM_MESH_NODAL_POLY_3D) {
    // Polyhedral mesh
    PDM_g_num_t ng_cell      = 0;
    PDM_g_num_t ng_face      = 0;
    PDM_g_num_t ng_vtx       = 0;
    int         dn_cell      = 0;
    int         dn_face      = 0;
    int         dn_edge      = 0;
    int         dn_vtx       = 0;
    int         n_face_group = 0;

    double      *dvtx_coord       = NULL;
    int         *dcell_face_idx   = NULL;
    PDM_g_num_t *dcell_face       = NULL;
    PDM_g_num_t *dface_cell       = NULL;
    int         *dface_vtx_idx    = NULL;
    PDM_g_num_t *dface_vtx        = NULL;
    int         *dface_group_idx  = NULL;
    PDM_g_num_t *dface_group      = NULL;

    PDM_poly_vol_gen(comm,
                     0.,
                     0.,
                     0.,
                     1.,
                     1.,
                     1.,
                     n_vtx_seg,
                     n_vtx_seg,
                     n_vtx_seg,
                     randomize,
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
                     &dface_vtx_idx,
                     &dface_vtx,
                     &dvtx_coord,
                     &dface_group_idx,
                     &dface_group);

    /* Generate dmesh */
    PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                          dn_cell,
                                          dn_face,
                                          dn_edge,
                                          dn_vtx,
                                          comm);

    PDM_dmesh_vtx_coord_set(dmesh,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);


    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_VTX,
                               dface_vtx,
                               dface_vtx_idx,
                               PDM_OWNERSHIP_KEEP);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               dface_cell,
                               NULL,
                               PDM_OWNERSHIP_KEEP);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_CELL_FACE,
                               dcell_face,
                               dcell_face_idx,
                               PDM_OWNERSHIP_KEEP);



    PDM_dmesh_bound_set(dmesh,
                        PDM_BOUND_TYPE_FACE,
                        n_face_group,
                        dface_group,
                        dface_group_idx,
                        PDM_OWNERSHIP_KEEP);

    PDM_dmesh_compute_distributions(dmesh);

    if (n_part > 0) {
      PDM_multipart_dmesh_set(*mpart, 0, dmesh);

      PDM_multipart_compute(*mpart);

      PDM_dmesh_free(dmesh);

      // Make edges great again (ugly AF)
      if (generate_edges) {
        int          *pn_face        = NULL;
        int          *pn_vtx         = NULL;
        int         **pface_vtx_idx  = NULL;
        int         **pface_vtx      = NULL;
        PDM_g_num_t **pface_ln_to_gn = NULL;
        PDM_g_num_t **pvtx_ln_to_gn  = NULL;
        PDM_malloc(pn_face       , n_part, int          );
        PDM_malloc(pn_vtx        , n_part, int          );
        PDM_malloc(pface_vtx_idx , n_part, int         *);
        PDM_malloc(pface_vtx     , n_part, int         *);
        PDM_malloc(pface_ln_to_gn, n_part, PDM_g_num_t *);
        PDM_malloc(pvtx_ln_to_gn , n_part, PDM_g_num_t *);

        int         **pface_edge_idx = NULL;
        int         **pface_edge     = NULL;
        int          *pn_edge        = NULL;
        int         **pedge_vtx      = NULL;
        PDM_g_num_t **pedge_ln_to_gn = NULL;

        for (int i_part = 0; i_part < n_part; i_part++) {
          PDM_multipart_part_connectivity_get(*mpart,
                                              0,
                                              i_part,
                                              PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                              &pface_vtx_idx[i_part],
                                              &pface_vtx    [i_part],
                                              PDM_OWNERSHIP_KEEP);

          pn_face[i_part] = PDM_multipart_part_ln_to_gn_get(*mpart,
                                                            0,
                                                            i_part,
                                                            PDM_MESH_ENTITY_FACE,
                                                            &pface_ln_to_gn[i_part],
                                                            PDM_OWNERSHIP_KEEP);

          pn_vtx[i_part] = PDM_multipart_part_ln_to_gn_get(*mpart,
                                                           0,
                                                           i_part,
                                                           PDM_MESH_ENTITY_VTX,
                                                           &pvtx_ln_to_gn[i_part],
                                                           PDM_OWNERSHIP_KEEP);
        }

        PDM_compute_face_edge_from_face_vtx(comm,
                                            n_part,
                                            pn_face,
                                            pn_vtx,
                                            pface_vtx_idx,
                                            pface_vtx,
                                            pface_ln_to_gn,
                                            pvtx_ln_to_gn,
                                            &pface_edge_idx,
                                            &pface_edge,
                                            &pn_edge,
                                            &pedge_vtx,
                                            &pedge_ln_to_gn);

        PDM_part_mesh_t *pm = (*mpart)->pmeshes[0].pmesh;
        for (int i_part = 0; i_part < n_part; i_part++) {
          PDM_part_mesh_connectivity_set(pm,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                         pface_edge    [i_part],
                                         pface_edge_idx[i_part],
                                         PDM_OWNERSHIP_KEEP);

          PDM_part_mesh_connectivity_set(pm,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                         pedge_vtx[i_part],
                                         NULL,
                                         PDM_OWNERSHIP_KEEP);

          PDM_part_mesh_n_entity_set(pm, i_part, PDM_MESH_ENTITY_EDGE, pn_edge[i_part]);
          PDM_part_mesh_entity_ln_to_gn_set(pm,
                                            i_part,
                                            PDM_MESH_ENTITY_EDGE,
                                            pedge_ln_to_gn[i_part],
                                            PDM_OWNERSHIP_KEEP);
        }
        PDM_free(pn_face       );
        PDM_free(pn_vtx        );
        PDM_free(pface_vtx_idx );
        PDM_free(pface_vtx     );
        PDM_free(pface_ln_to_gn);
        PDM_free(pvtx_ln_to_gn );
        PDM_free(pface_edge_idx);
        PDM_free(pface_edge    );
        PDM_free(pn_edge       );
        PDM_free(pedge_vtx     );
        PDM_free(pedge_ln_to_gn);
      }
    }
    else {
      *out_dmesh = dmesh;
    }
  }
  else {
    PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                          n_vtx_seg,
                                                          n_vtx_seg,
                                                          n_vtx_seg,
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
                                       PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                       PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

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
      else {
        int *face_vtx_idx = NULL;
        int *face_vtx     = NULL;
        PDM_multipart_part_connectivity_get(*mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                            &face_vtx_idx,
                                            &face_vtx,
                                            PDM_OWNERSHIP_KEEP);
        assert(face_vtx != NULL);
        PDM_part_mesh_connectivity_set(pmesh,
                                       i_part,
                                       PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                       face_vtx,
                                       face_vtx_idx,
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


static void
_compute_itp_field
(
 int     n_vtx,
 double *vtx_coord,
 double *vtx_field
 )
{
  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    vtx_field[i_vtx] = cos(5*(vtx_coord[3*i_vtx+1] + vtx_coord[3*i_vtx+2]));
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
  PDM_Mesh_nodal_elt_t  elt_type       = PDM_MESH_NODAL_HEXA8;
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


  /**
   * Compute scalar field
   */
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


  /**
   * Create Isosurface instance
   */
  PDM_isosurface_t *isos = PDM_isosurface_create(comm, 3);

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
  if (n_part > 0) {
    // Partitioned
    if (use_part_mesh) {
      PDM_isosurface_part_mesh_set(isos, pmesh);
    }
    else {
      PDM_isosurface_n_part_set(isos, n_part);

      for (int i_part = 0; i_part < n_part; i_part++) {

        // Connectivities
        int *cell_face_idx = NULL;
        int *cell_face     = NULL;
        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                            &cell_face_idx,
                                            &cell_face,
                                            PDM_OWNERSHIP_KEEP);

        int *face_edge_idx = NULL;
        int *face_edge     = NULL;
        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                            &face_edge_idx,
                                            &face_edge,
                                            PDM_OWNERSHIP_KEEP);

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
                                        PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                        cell_face_idx,
                                        cell_face);

        PDM_isosurface_connectivity_set(isos,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                        face_edge_idx,
                                        face_edge);

        PDM_isosurface_connectivity_set(isos,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        NULL,
                                        edge_vtx);

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
        PDM_g_num_t *cell_ln_to_gn = NULL;
        int n_cell = PDM_multipart_part_ln_to_gn_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_MESH_ENTITY_CELL,
                                                     &cell_ln_to_gn,
                                                     PDM_OWNERSHIP_KEEP);

        PDM_g_num_t *face_ln_to_gn = NULL;
        int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_MESH_ENTITY_FACE,
                                                     &face_ln_to_gn,
                                                     PDM_OWNERSHIP_KEEP);

        PDM_g_num_t *edge_ln_to_gn = NULL;
        int n_edge = PDM_multipart_part_ln_to_gn_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_MESH_ENTITY_EDGE,
                                                     &edge_ln_to_gn,
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
                                    PDM_MESH_ENTITY_CELL,
                                    n_cell,
                                    cell_ln_to_gn);

        PDM_isosurface_ln_to_gn_set(isos,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    n_face,
                                    face_ln_to_gn);

        PDM_isosurface_ln_to_gn_set(isos,
                                    i_part,
                                    PDM_MESH_ENTITY_EDGE,
                                    n_edge,
                                    edge_ln_to_gn);

        PDM_isosurface_ln_to_gn_set(isos,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    n_vtx,
                                    vtx_ln_to_gn);

        // Groups
        int          n_surface             = 0;
        int         *surface_face_idx      = NULL;
        int         *surface_face          = NULL;
        PDM_g_num_t *surface_face_ln_to_gn = NULL;
        PDM_multipart_group_get(mpart,
                                0,
                                i_part,
                                PDM_MESH_ENTITY_FACE,
                                &n_surface,
                                &surface_face_idx,
                                &surface_face,
                                &surface_face_ln_to_gn,
                                PDM_OWNERSHIP_KEEP);

        PDM_isosurface_group_set(isos,
                                 i_part,
                                 PDM_MESH_ENTITY_FACE,
                                 n_surface,
                                 surface_face_idx,
                                 surface_face,
                                 surface_face_ln_to_gn);
      }
    }
  }
  else {
    // Block-distributed
    if (use_part_mesh) {
      PDM_isosurface_dmesh_set(isos, dmesh);
    }
    else {
      for (int i_entity = PDM_MESH_ENTITY_CELL; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
        PDM_g_num_t *distrib = NULL;
        PDM_dmesh_distrib_get     (dmesh, i_entity, &distrib);
        PDM_isosurface_distrib_set(isos,  i_entity,  distrib);
      }

      int         *dcell_face_idx = NULL;
      PDM_g_num_t *dcell_face     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                 &dcell_face,
                                 &dcell_face_idx,
                                 PDM_OWNERSHIP_KEEP);

      PDM_isosurface_dconnectivity_set(isos,
                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                       dcell_face_idx,
                                       dcell_face);

      int         *dface_edge_idx = NULL;
      PDM_g_num_t *dface_edge     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                 &dface_edge,
                                 &dface_edge_idx,
                                 PDM_OWNERSHIP_KEEP);
      PDM_isosurface_dconnectivity_set(isos,
                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                       dface_edge_idx,
                                       dface_edge);

      int         *dface_vtx_idx = NULL;
      PDM_g_num_t *dface_vtx     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                 &dface_vtx,
                                 &dface_vtx_idx,
                                 PDM_OWNERSHIP_KEEP);
      PDM_isosurface_dconnectivity_set(isos,
                                       PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                       dface_vtx_idx,
                                       dface_vtx);

      int         *dedge_vtx_idx = NULL;
      PDM_g_num_t *dedge_vtx     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                 &dedge_vtx,
                                 &dedge_vtx_idx,
                                 PDM_OWNERSHIP_KEEP);

      PDM_isosurface_dconnectivity_set(isos,
                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                       NULL,
                                       dedge_vtx);

      double *dvtx_coord = NULL;
      PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_KEEP);
      PDM_isosurface_dvtx_coord_set(isos, dvtx_coord);

      int         *dsurface_face_idx = NULL;
      PDM_g_num_t *dsurface_face     = NULL;
      int n_surface = PDM_dmesh_bound_get(dmesh,
                                          PDM_BOUND_TYPE_FACE,
                                          &dsurface_face,
                                          &dsurface_face_idx,
                                          PDM_OWNERSHIP_KEEP);
      PDM_isosurface_dgroup_set(isos,
                                PDM_MESH_ENTITY_FACE,
                                n_surface,
                                dsurface_face_idx,
                                dsurface_face);
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

  int n_iso = iso2 + 1;

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
  for (int i_iso = 0; i_iso < n_iso; i_iso++) {
    PDM_isosurface_enable_part_to_part(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_VTX,
                                       0);

    PDM_isosurface_enable_part_to_part(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_EDGE,
                                       0);

    PDM_isosurface_enable_part_to_part(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_FACE,
                                       0);
  }

  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_reset  (isos, iso1);
  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_compute(isos, iso2);
  PDM_isosurface_compute(isos, -1);
  // PDM_isosurface_reset  (isos, -1);
  // PDM_isosurface_compute(isos, -1);


  /*
   *  Interpolate field
   */
  double **itp_field  = NULL;
  // double  *itp_dfield = NULL;

  double ***iso_itp_field   = NULL;
  double ***iso_parent_cell = NULL;
  double ***iso_parent_face = NULL;
  // double  **iso_itp_dfield = NULL;

  if (n_part > 0) {
    // Partitioned
    PDM_malloc(iso_itp_field,   n_iso, double **);
    PDM_malloc(iso_parent_cell, n_iso, double **);
    PDM_malloc(iso_parent_face, n_iso, double **);

    PDM_malloc(itp_field, n_part, double *);
    for (int i_part = 0; i_part < n_part; i_part++) {

      double *vtx_coord = NULL;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx_coord,
                                                   PDM_OWNERSHIP_KEEP);
      PDM_malloc(itp_field[i_part], n_vtx, double);

      _compute_itp_field(n_vtx, vtx_coord, itp_field[i_part]);
    }


    for (int i_iso = 0; i_iso < n_iso; i_iso++) {

      PDM_malloc(iso_itp_field[i_iso], n_part, double *);

      if (local) {
        // Local
        for (int i_part = 0; i_part < n_part; i_part++) {
          int    *iso_vtx_parent_idx;
          double *iso_vtx_parent_weight;
          int iso_n_vtx = PDM_isosurface_vtx_parent_weight_get(isos,
                                                               i_iso,
                                                               i_part,
                                                               &iso_vtx_parent_idx,
                                                               &iso_vtx_parent_weight,
                                                               PDM_OWNERSHIP_KEEP);

          int *iso_vtx_parent;
          PDM_isosurface_local_parent_get(isos,
                                          i_iso,
                                          i_part,
                                          PDM_MESH_ENTITY_VTX,
                                          &iso_vtx_parent_idx,
                                          &iso_vtx_parent,
                                          PDM_OWNERSHIP_KEEP);

          PDM_malloc(iso_itp_field[i_iso][i_part], iso_n_vtx, double);
          for (int i_vtx = 0; i_vtx < iso_n_vtx; i_vtx++) {
            iso_itp_field[i_iso][i_part][i_vtx] = 0.;
            for (int i = iso_vtx_parent_idx[i_vtx]; i < iso_vtx_parent_idx[i_vtx+1]; i++) {
              int i_parent = iso_vtx_parent[i_vtx] - 1;
              iso_itp_field[i_iso][i_part][i_vtx] += iso_vtx_parent_weight[i] * itp_field[i_part][i_parent];
            }
          }
        } // End loop on parts
      } // End if LOCAL

      else {
        // Reequilibrate

        // Interpolate a vtx-based field
        PDM_part_to_part_t *ptp_vtx = NULL;
        PDM_isosurface_part_to_part_get(isos,
                                        i_iso,
                                        PDM_MESH_ENTITY_VTX,
                                        &ptp_vtx,
                                        PDM_OWNERSHIP_KEEP);

        double **recv_vtx_field = NULL;
        int request_vtx = -1;
        PDM_part_to_part_reverse_iexch(ptp_vtx,
                                       PDM_MPI_COMM_KIND_P2P,
                                       PDM_STRIDE_CST_INTERLACED,
                                       PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                       1,
                                       sizeof(double),
                                       NULL,
                      (const void  **) itp_field,
                                       NULL,
                      (      void ***) &recv_vtx_field,
                                       &request_vtx);

        PDM_part_to_part_reverse_iexch_wait(ptp_vtx, request_vtx);

        for (int i_part = 0; i_part < n_part; i_part++) {
          int    *iso_vtx_parent_idx;
          double *iso_vtx_parent_weight;
          int iso_n_vtx = PDM_isosurface_vtx_parent_weight_get(isos,
                                                               i_iso,
                                                               i_part,
                                                               &iso_vtx_parent_idx,
                                                               &iso_vtx_parent_weight,
                                                               PDM_OWNERSHIP_KEEP);

          PDM_malloc(iso_itp_field[i_iso][i_part], iso_n_vtx, double);
          for (int i_vtx = 0; i_vtx < iso_n_vtx; i_vtx++) {
            iso_itp_field[i_iso][i_part][i_vtx] = 0.;
            for (int i = iso_vtx_parent_idx[i_vtx]; i < iso_vtx_parent_idx[i_vtx+1]; i++) {
              iso_itp_field[i_iso][i_part][i_vtx] += iso_vtx_parent_weight[i] * recv_vtx_field[i_part][i];
            }
          }

          PDM_free(recv_vtx_field[i_part]);
        } // End loop on parts
        PDM_free(recv_vtx_field);


        // Exchange parent cell gnums just for fun
        PDM_part_to_part_t *ptp_face = NULL;
        PDM_isosurface_part_to_part_get(isos,
                                        i_iso,
                                        PDM_MESH_ENTITY_FACE,
                                        &ptp_face,
                                        PDM_OWNERSHIP_KEEP);

        PDM_g_num_t **pcell_ln_to_gn;
        PDM_malloc(pcell_ln_to_gn, n_part, PDM_g_num_t *);
        for (int i_part = 0; i_part < n_part; i_part++) {
          PDM_multipart_part_ln_to_gn_get(mpart,
                                          0,
                                          i_part,
                                          PDM_MESH_ENTITY_CELL,
                                          &pcell_ln_to_gn[i_part],
                                          PDM_OWNERSHIP_KEEP);
        }

        PDM_g_num_t **recv_face_field = NULL;
        int request_face = -1;
        PDM_part_to_part_reverse_iexch(ptp_face,
                                       PDM_MPI_COMM_KIND_P2P,
                                       PDM_STRIDE_CST_INTERLACED,
                                       PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                       1,
                                       sizeof(PDM_g_num_t),
                                       NULL,
                      (const void  **) pcell_ln_to_gn,
                                       NULL,
                      (      void ***) &recv_face_field,
                                       &request_face);

        PDM_part_to_part_reverse_iexch_wait(ptp_face, request_face);
        PDM_free(pcell_ln_to_gn);

        // convert gnums to doubles for vtk export
        PDM_malloc(iso_parent_cell[i_iso], n_part, double *);
        for (int i_part = 0; i_part < n_part; i_part++) {
          PDM_g_num_t *iso_face_ln_to_gn;
          int iso_n_face = PDM_isosurface_ln_to_gn_get(isos,
                                                       i_iso,
                                                       i_part,
                                                       PDM_MESH_ENTITY_FACE,
                                                       &iso_face_ln_to_gn,
                                                       PDM_OWNERSHIP_KEEP);

          PDM_malloc(iso_parent_cell[i_iso][i_part], iso_n_face, double);
          for (int i = 0; i < iso_n_face; i++) {
            iso_parent_cell[i_iso][i_part][i] = recv_face_field[i_part][i];
          }
          PDM_free(recv_face_field[i_part]);
        }
        PDM_free(recv_face_field);



        // Exchange parent face gnums just for fun
        PDM_part_to_part_t *ptp_edge = NULL;
        PDM_isosurface_part_to_part_get(isos,
                                        i_iso,
                                        PDM_MESH_ENTITY_EDGE,
                                        &ptp_edge,
                                        PDM_OWNERSHIP_KEEP);

        PDM_g_num_t **pface_ln_to_gn;
        PDM_malloc(pface_ln_to_gn, n_part, PDM_g_num_t *);
        for (int i_part = 0; i_part < n_part; i_part++) {
          PDM_multipart_part_ln_to_gn_get(mpart,
                                          0,
                                          i_part,
                                          PDM_MESH_ENTITY_FACE,
                                          &pface_ln_to_gn[i_part],
                                          PDM_OWNERSHIP_KEEP);
        }

        PDM_g_num_t **recv_edge_field = NULL;
        int request_edge = -1;
        PDM_part_to_part_reverse_iexch(ptp_edge,
                                       PDM_MPI_COMM_KIND_P2P,
                                       PDM_STRIDE_CST_INTERLACED,
                                       PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                       1,
                                       sizeof(PDM_g_num_t),
                                       NULL,
                      (const void  **) pface_ln_to_gn,
                                       NULL,
                      (      void ***) &recv_edge_field,
                                       &request_edge);

        PDM_part_to_part_reverse_iexch_wait(ptp_edge, request_edge);
        PDM_free(pface_ln_to_gn);

        // convert gnums to doubles for vtk export
        PDM_malloc(iso_parent_face[i_iso], n_part, double *);
        for (int i_part = 0; i_part < n_part; i_part++) {
          PDM_g_num_t *iso_edge_ln_to_gn;
          int iso_n_edge = PDM_isosurface_ln_to_gn_get(isos,
                                                       i_iso,
                                                       i_part,
                                                       PDM_MESH_ENTITY_EDGE,
                                                       &iso_edge_ln_to_gn,
                                                       PDM_OWNERSHIP_KEEP);

          PDM_malloc(iso_parent_face[i_iso][i_part], iso_n_edge, double);
          for (int i = 0; i < iso_n_edge; i++) {
            iso_parent_face[i_iso][i_part][i] = recv_edge_field[i_part][i];
          }
          PDM_free(recv_edge_field[i_part]);
        }
        PDM_free(recv_edge_field);

      } // End if REEQUILIBRATE


    }


  }
  else {
    // Block-distributed
    // TODO...
    printf("TODO: field interplation in block-distributed\n");
  }


  /*
   *  Visu isosurfaces
   */
  if (visu) {
    char out_name[999];
    for (int i_iso = 0; i_iso < n_iso; i_iso++) {

      if (n_part > 0) {
        // Partitioned
        for (int i_part = 0; i_part < n_part; i_part++) {
          double *iso_vtx_coord = NULL;
          int iso_n_vtx = PDM_isosurface_vtx_coord_get(isos,
                                                       i_iso,
                                                       i_part,
                                                       &iso_vtx_coord,
                                                       PDM_OWNERSHIP_KEEP);

          int *iso_face_vtx_idx = NULL;
          int *iso_face_vtx     = NULL;
          int iso_n_face = PDM_isosurface_connectivity_get(isos,
                                                           i_iso,
                                                           i_part,
                                                           PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                           &iso_face_vtx_idx,
                                                           &iso_face_vtx,
                                                           PDM_OWNERSHIP_KEEP);

          int *iso_edge_vtx_idx = NULL;
          int *iso_edge_vtx     = NULL;
          int iso_n_edge = PDM_isosurface_connectivity_get(isos,
                                                           i_iso,
                                                           i_part,
                                                           PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                           &iso_edge_vtx_idx,
                                                           &iso_edge_vtx,
                                                           PDM_OWNERSHIP_KEEP);

          PDM_g_num_t *iso_vtx_ln_to_gn = NULL;
          PDM_isosurface_ln_to_gn_get(isos,
                                      i_iso,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &iso_vtx_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

          PDM_g_num_t *iso_face_ln_to_gn = NULL;
          PDM_isosurface_ln_to_gn_get(isos,
                                      i_iso,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &iso_face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

          PDM_g_num_t *iso_edge_ln_to_gn = NULL;
          PDM_isosurface_ln_to_gn_get(isos,
                                      i_iso,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &iso_edge_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

          int *isovalue_face_idx = NULL;
          int _n_isovalues = PDM_isosurface_isovalue_entity_idx_get(isos,
                                                                    i_iso,
                                                                    i_part,
                                                                    PDM_MESH_ENTITY_FACE,
                                                                    &isovalue_face_idx,
                                                                    PDM_OWNERSHIP_KEEP);

          double *iso_face_isovalue = NULL;
          PDM_malloc(iso_face_isovalue, iso_n_face, double);
          for (int i_isovalue = 0; i_isovalue < _n_isovalues; i_isovalue++) {
            for (int i_face = isovalue_face_idx[i_isovalue]; i_face < isovalue_face_idx[i_isovalue+1]; i_face++) {
              iso_face_isovalue[i_face] = i_isovalue;
            }
          }


          sprintf(out_name, "isosurface_3d_ngon_iso_%d_part_%d.vtk", i_iso, i_rank*n_part + i_part);
          PDM_vtk_write_polydata_field(out_name,
                                       iso_n_vtx,
                                       iso_vtx_coord,
                                       iso_vtx_ln_to_gn,
                                       iso_n_face,
                                       iso_face_vtx_idx,
                                       iso_face_vtx,
                                       iso_face_ln_to_gn,
                                       "parent_cell",//"isovalue",
                                       iso_parent_cell[i_iso][i_part],//iso_face_isovalue,
                                       "itp_field",
                                       iso_itp_field[i_iso][i_part]);
          PDM_free(iso_face_isovalue);

          int *isovalue_edge_idx = NULL;
          _n_isovalues = PDM_isosurface_isovalue_entity_idx_get(isos,
                                                                i_iso,
                                                                i_part,
                                                                PDM_MESH_ENTITY_EDGE,
                                                                &isovalue_edge_idx,
                                                                PDM_OWNERSHIP_KEEP);

          double *iso_edge_isovalue = NULL;
          double *iso_edge_group    = NULL;
          PDM_malloc(iso_edge_isovalue, iso_n_edge, double);
          PDM_malloc(iso_edge_group   , iso_n_edge, double);
          for (int i_isovalue = 0; i_isovalue < _n_isovalues; i_isovalue++) {
            for (int i_edge = isovalue_edge_idx[i_isovalue]; i_edge < isovalue_edge_idx[i_isovalue+1]; i_edge++) {
              iso_edge_isovalue[i_edge] = i_isovalue;
            }
          }


          int          n_group             = 0;
          int         *group_edge_idx      = NULL;
          int         *group_edge          = NULL;
          PDM_g_num_t *group_edge_ln_to_gn = NULL;
          n_group = PDM_isosurface_group_get(isos,
                                             i_iso,
                                             i_part,
                                             PDM_MESH_ENTITY_EDGE,
                                             &group_edge_idx,
                                             &group_edge,
                                             &group_edge_ln_to_gn,
                                             PDM_OWNERSHIP_KEEP);
          for (int i_group = 0; i_group < n_group; i_group++) {
            for (int i = group_edge_idx[i_group]; i < group_edge_idx[i_group+1]; i++) {
              iso_edge_group[group_edge[i]-1] = i_group + 1;
            }
          }

          sprintf(out_name, "isosurface_3d_ngon_iso_edge_%d_part_%d.vtk", i_iso, i_rank*n_part + i_part);
          const char   *elt_field_name [] = {"i_isovalue", "i_group", "parent_face"};
          const double *elt_field_value[] = {iso_edge_isovalue, iso_edge_group, iso_parent_face[i_iso][i_part]};

          const char   *vtx_field_name [] = {"itp_field"};
          const double *vtx_field_value[] = {(const double *) iso_itp_field[i_iso][i_part]};

          PDM_vtk_write_std_elements_ho_with_vtx_field(out_name,
                                                       1,
                                                       iso_n_vtx,
                                                       iso_vtx_coord,
                                                       iso_vtx_ln_to_gn,
                                                       PDM_MESH_NODAL_BAR2,
                                                       iso_n_edge,
                                                       iso_edge_vtx,
                                                       iso_edge_ln_to_gn,
                                                       3,
                                                       elt_field_name,
                                                       elt_field_value,
                                                       1,
                                                       vtx_field_name,
                                                       vtx_field_value);
          PDM_free(iso_edge_isovalue);
          PDM_free(iso_edge_group   );
        }
      }
      else {
        // > Vertices
        double *dvtx_coords = NULL;
        PDM_g_num_t dn_vtx = PDM_isosurface_dvtx_coord_get(isos, i_iso, &dvtx_coords, PDM_OWNERSHIP_KEEP);
        if (1) {
          log_trace("dn_vtx = "PDM_FMT_G_NUM"\n", dn_vtx);
          PDM_log_trace_array_double(dvtx_coords, 3*dn_vtx, "dvtx_coords ::");
        }

        // > Edges
        int         *dedge_vtx_idx = NULL;
        PDM_g_num_t *dedge_vtx     = NULL;
        PDM_g_num_t dn_edge = PDM_isosurface_dconnectivity_get(isos, i_iso, 
                                                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                              &dedge_vtx_idx,
                                                              &dedge_vtx,
                                                               PDM_OWNERSHIP_KEEP);
        if (1) {
          log_trace("dn_edge = "PDM_FMT_G_NUM"\n", dn_edge);
          PDM_log_trace_array_long(dedge_vtx, 2*dn_edge, "dedge_vtx ::");
        }

        // > Faces
        int         *dface_vtx_idx = NULL;
        PDM_g_num_t *dface_vtx     = NULL;
        PDM_g_num_t dn_face = PDM_isosurface_dconnectivity_get(isos, i_iso, 
                                                               PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                              &dface_vtx_idx,
                                                              &dface_vtx,
                                                               PDM_OWNERSHIP_KEEP);
        if (1) {
          log_trace("dn_face = "PDM_FMT_G_NUM"\n", dn_face);
          int size_face_vtx = dface_vtx_idx[dn_face];
          PDM_log_trace_array_int (dface_vtx_idx, dn_face      , "dface_vtx_idx ::");
          PDM_log_trace_array_long(dface_vtx    , size_face_vtx, "dface_vtx     ::");
        }

        // > Groups
        int         *dedge_group_idx = NULL;
        PDM_g_num_t *dedge_group     = NULL;
        int n_group = PDM_isosurface_dgroup_get(isos, i_iso, 
                                                PDM_MESH_ENTITY_EDGE,
                                               &dedge_group_idx,
                                               &dedge_group,
                                                PDM_OWNERSHIP_KEEP);
        if (1) {
          log_trace("n_group = %d\n", n_group);
          int size_edge_group = dedge_group_idx[n_group];
          PDM_log_trace_array_int (dedge_group_idx, n_group        , "dedge_group_idx ::");
          PDM_log_trace_array_long(dedge_group    , size_edge_group, "dedge_group ::");
        }

        /**
         * Prepare dmesh_nodal for ouput
         */
        PDM_g_num_t  n_entity[3] = {dn_vtx, dn_edge, dn_face};
        PDM_g_num_t gn_entity[3];
        PDM_MPI_Allreduce(n_entity, gn_entity, 3, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
        PDM_dmesh_nodal_t *iso_dmn = PDM_DMesh_nodal_create(comm,
                                                            2,
                                                            gn_entity[0],
                                                            0,
                                                            gn_entity[2],
                                                            gn_entity[1]);

        PDM_DMesh_nodal_coord_set(iso_dmn, dn_vtx, dvtx_coords, PDM_OWNERSHIP_USER);

        int edge_section = PDM_DMesh_nodal_section_add(iso_dmn,
                                                       PDM_GEOMETRY_KIND_RIDGE,
                                                       PDM_MESH_NODAL_BAR2);
        PDM_DMesh_nodal_section_std_set(iso_dmn,
                                        PDM_GEOMETRY_KIND_RIDGE,
                                        edge_section,
                                        dn_edge,
                                        dedge_vtx,
                                        PDM_OWNERSHIP_USER);

        int face_section = PDM_DMesh_nodal_section_add(iso_dmn,
                                                       PDM_GEOMETRY_KIND_SURFACIC,
                                                       PDM_MESH_NODAL_POLY_2D);

        PDM_DMesh_nodal_section_poly2d_set(iso_dmn,
                                           PDM_GEOMETRY_KIND_SURFACIC,
                                           face_section,
                                           dn_face,
                                           dface_vtx_idx,
                                           dface_vtx,
                                           PDM_OWNERSHIP_USER);

        PDM_DMesh_nodal_section_group_elmt_set(iso_dmn,
                                               PDM_GEOMETRY_KIND_RIDGE,
                                               n_group,
                                               dedge_group_idx,
                                               dedge_group,
                                               PDM_OWNERSHIP_USER);

        sprintf(out_name, "isosurface_3d_ngon_iso_edge_%d.vtk", i_iso);
        PDM_dmesh_nodal_dump_vtk(iso_dmn, PDM_GEOMETRY_KIND_RIDGE, out_name);
        sprintf(out_name, "isosurface_3d_ngon_iso_%d.vtk", i_iso);
        PDM_dmesh_nodal_dump_vtk(iso_dmn, PDM_GEOMETRY_KIND_SURFACIC, out_name);

        PDM_DMesh_nodal_free(iso_dmn);
      }

    }
  }


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
      PDM_free(itp_field [i_part]);
    }
    PDM_free(piso_field);
    PDM_free(itp_field );

    for (int i_iso = 0; i_iso < n_iso; i_iso++) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        PDM_free(iso_itp_field  [i_iso][i_part]);
        PDM_free(iso_parent_cell[i_iso][i_part]);
        PDM_free(iso_parent_face[i_iso][i_part]);
      }
      PDM_free(iso_itp_field  [i_iso]);
      PDM_free(iso_parent_cell[i_iso]);
      PDM_free(iso_parent_face[i_iso]);
    }
    PDM_free(iso_itp_field  );
    PDM_free(iso_parent_cell);
    PDM_free(iso_parent_face);
  }
  else {
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
