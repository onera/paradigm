#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static PDM_multipart_t *
_gen_mesh
(
 PDM_MPI_Comm          comm,
 const char           *filename,
 int                   n_part,
 PDM_g_num_t           n_vtx_seg,
 int                   randomize,
 PDM_Mesh_nodal_elt_t  elt_type
 )
{
  int n_domain = 1;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                               &n_part,
                                                PDM_FALSE,
                                                PDM_SPLIT_DUAL_WITH_HILBERT,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  if (filename != NULL) {
    PDM_dmesh_nodal_t *dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                                          filename,
                                                          0,
                                                          0);

    if (0) {
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC,  "isosurface_3d_ngon_volume");
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "isosurface_3d_ngon_surface");
    }


    PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
    PDM_multipart_compute(mpart);

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
    free(dcell_face_idx);
    free(dcell_face    );

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

    PDM_dmesh_bound_set(dmesh,
                        PDM_BOUND_TYPE_FACE,
                        n_face_group,
                        dface_group,
                        dface_group_idx,
                        PDM_OWNERSHIP_KEEP);

    PDM_multipart_dmesh_set(mpart, 0, dmesh);

    PDM_multipart_compute(mpart);

    PDM_dmesh_free(dmesh);
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
                                                          PDM_OWNERSHIP_KEEP);

    PDM_dcube_nodal_gen_random_factor_set(dcube, (double) randomize);

    PDM_dcube_nodal_gen_build(dcube);

    PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

    PDM_dmesh_nodal_generate_distribution(dmn);

    PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
    PDM_multipart_compute(mpart);

    PDM_dcube_nodal_gen_free(dcube);
  }

  return mpart;
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

// static void
// _check_result
// (
// )


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
     "  -h                           This message.\n\n"
     "  -in          <filename>      Mesh file name (Gamma Mesh Format).\n\n"
     "  -sol         <filename>      Solution file name (Gamma Mesh Format).\n\n"
     "  -visu                        Enable exports for visualization.\n\n"
     "  -n_isovalues <n>             Number of isovalues.\n\n"
     "  -isovalues   <v1 v2 ... vn>  Isovalues.\n\n"
     "  -elt_type    <t>             Volume element type (only for automatically generated mesh).\n\n"
     "  -randomize                   Enable mesh randomization (only for automatically generated mesh).\n\n"
     "  -n           <n>             Number of vertices per side (only for automatically generated mesh).\n\n"
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
  PDM_g_num_t           *n_vtx_seg
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
      *isovalues = malloc(sizeof(double) * (*n_isovalues));
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

/**
 * TODO:
 *  * block-ditributed mesh
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
  char                 *mesh_name   = NULL;
  char                 *sol_name    = NULL;
  int                   n_part      = 1;
  int                   visu        = 0;
  int                   n_isovalues = 1;
  double               *isovalues   = NULL;
  PDM_Mesh_nodal_elt_t  elt_type    = PDM_MESH_NODAL_HEXA8;
  PDM_g_num_t           n_vtx_seg   = 10;
  int                   randomize   = 0;
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
             &n_vtx_seg);

  if (n_part < 1) {
    PDM_error(__FILE__, __LINE__, 0, "Block-distributed not yet available.");
  }

  if (isovalues == NULL) {
    n_isovalues = 1;
    isovalues = malloc(sizeof(double) * n_isovalues);
    isovalues[0] = 0.;
  }

  // Generate/read mesh
  PDM_multipart_t *mpart = _gen_mesh(comm,
                                     mesh_name,
                                     n_part,
                                     n_vtx_seg,
                                     randomize,
                                     elt_type);

  // Compute scalar field
  double **iso_field = NULL;
  if (n_part > 0) {
    // Partitioned
    iso_field = malloc(sizeof(double *) * n_part);
    for (int i_part = 0; i_part < n_part; i_part++) {

      double *vtx_coord = NULL;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx_coord,
                                                   PDM_OWNERSHIP_KEEP);

      iso_field[i_part] = malloc(sizeof(double) * n_vtx);
      _compute_iso_field(n_vtx, vtx_coord, iso_field[i_part]);
    }
  }
  else {
    // Block-distributed
    // TODO
  }


  // Create Isosurface instance
  PDM_isosurface_t *isos = PDM_isosurface_create(comm, 3);

  // TODO:
  // - PDM_isosurface_redistribution_set


  /* Set mesh */
  if (n_part > 0) {
    // Partitioned
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
  else {
    // Block-distributed
    // TODO
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
  // double field_isovalues[1] = {0.1};
  int iso2 = PDM_isosurface_add(isos,
                                PDM_ISO_SURFACE_KIND_FIELD,
                                n_isovalues,// 1,
                                isovalues);// field_isovalues);

  if (n_part > 0) {
    // Partitioned
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_isosurface_field_set(isos,
                               iso2,
                               i_part,
                               iso_field[i_part]);
    }
  }
  else {
    // Block-distributed
    // TODO
  }



  /*
   *  Compute isosurface
   */
  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_reset  (isos, iso1);
  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_compute(isos, iso2);


  /*
   *  Visu isosurfaces
   */
  int n_iso = iso2 + 1;
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
                                                                    &isovalue_face_idx);

          int *iso_face_isovalue = malloc(sizeof(int) * iso_n_face);
          for (int i_isovalue = 0; i_isovalue < _n_isovalues; i_isovalue++) {
            for (int i_face = isovalue_face_idx[i_isovalue]; i_face < isovalue_face_idx[i_isovalue+1]; i_face++) {
              iso_face_isovalue[i_face] = i_isovalue;
            }
          }


          sprintf(out_name, "isosurface_3d_ngon_iso_%d_part_%d.vtk", i_iso, i_rank*n_part + i_part);
          PDM_vtk_write_polydata(out_name,
                                 iso_n_vtx,
                                 iso_vtx_coord,
                                 iso_vtx_ln_to_gn,
                                 iso_n_face,
                                 iso_face_vtx_idx,
                                 iso_face_vtx,
                                 NULL,//iso_face_ln_to_gn,
                                 iso_face_isovalue);
          free(iso_face_isovalue);

          // TODO: edges
          int *isovalue_edge_idx = NULL;
          _n_isovalues = PDM_isosurface_isovalue_entity_idx_get(isos,
                                                                i_iso,
                                                                i_part,
                                                                PDM_MESH_ENTITY_EDGE,
                                                                &isovalue_edge_idx);

          int *iso_edge_isovalue = malloc(sizeof(int) * iso_n_edge);
          int *iso_edge_group    = malloc(sizeof(int) * iso_n_edge);
          for (int i_isovalue = 0; i_isovalue < _n_isovalues; i_isovalue++) {
            for (int i_edge = isovalue_edge_idx[i_isovalue]; i_edge < isovalue_edge_idx[i_isovalue+1]; i_edge++) {
              iso_edge_isovalue[i_edge] = i_isovalue;
            }
          }


          int          n_group             = 0;
          int         *group_edge_idx      = NULL;
          int         *group_edge          = NULL;
          PDM_g_num_t *group_edge_ln_to_gn = NULL;
          PDM_isosurface_group_get(isos,
                                   i_iso,
                                   i_part,
                                   PDM_MESH_ENTITY_EDGE,
                                   &n_group,
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
          const char *field_name [] = {"i_isovalue", "i_group"};
          const int  *field_value[] = {iso_edge_isovalue, iso_edge_group};

          PDM_vtk_write_std_elements(out_name,
                                     iso_n_vtx,
                                     iso_vtx_coord,
                                     iso_vtx_ln_to_gn,
                                     PDM_MESH_NODAL_BAR2,
                                     iso_n_edge,
                                     iso_edge_vtx,
                                     iso_edge_ln_to_gn,
                                     2,
                                     field_name,
                                     field_value);
          free(iso_edge_isovalue);
          free(iso_edge_group   );
        }
      }
      else {
        // Block-distributed
        // TODO
      }


    }
  }


  /* Free memory */
  PDM_isosurface_free(isos);

  PDM_multipart_free(mpart);

  free(isovalues);

  if (n_part > 0) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      free(iso_field[i_part]);
    }
    free(iso_field);
  }


  if (i_rank == 0) {
    printf("End :D\n");
    fflush(stdout);
  }



  PDM_MPI_Finalize();


  return EXIT_SUCCESS;
}
