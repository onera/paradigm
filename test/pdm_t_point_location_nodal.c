#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_generate_mesh.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_point_location.h"
#include "pdm_printf.h"
#include "pdm_priv.h"
#include "pdm_vtk.h"
#include "pdm.h"

/*============================================================================
 * Private function definitions
 *============================================================================*/

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n         <int>    Number of vertices on the cube side.\n\n"
     "  -l         <float>  Cube length.\n\n"
     "  -n_part    <int>    Number of partitions par process.\n\n"
     "  -post               Enable outputs for visualization.\n\n"
     "  -parmetis           Call ParMETIS.\n\n"
     "  -pt-scocth          Call PT-Scotch.\n\n"
     "  -t         <int>    Type of mesh elements.\n\n"
     "  -f         <str>    Path to mesh file.\n\n"
     "  -h                  This message.\n\n");

  exit(exit_code);
}

static void
_read_args
(
  int                    argc,
  char                 **argv,
  PDM_g_num_t           *n_vtx_seg,
  double                *length,
  int                   *n_part,
  int                   *post,
  PDM_split_dual_t      *part_method,
  PDM_Mesh_nodal_elt_t  *elt_type,
  char                 **filename
)
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
      *part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_generate_mesh
(
        PDM_MPI_Comm            comm,
        int                     n_part,
        PDM_split_dual_t        part_method,
        PDM_Mesh_nodal_elt_t    elt_type,
        PDM_g_num_t             n_vtx_seg,
  const char                   *filename,
        PDM_part_mesh_nodal_t **out_pmn
)
{
  if (filename != NULL) {
    // Read mesh from file
    *out_pmn = PDM_generate_mesh_nodal_from_file(comm,
                                                 n_part,
                                                 part_method,
                                                 filename);
  }

  else {
    // Generate procedural mesh

    int dim = PDM_Mesh_nodal_elt_dim_get(elt_type);

    if (dim == 2) {
      *out_pmn = PDM_generate_mesh_rectangle(comm,
                                             elt_type,
                                             1,
                                             NULL,
                                             -1,
                                             -1,
                                             0,
                                             2,
                                             2,
                                             n_vtx_seg,
                                             n_vtx_seg,
                                             n_part,
                                             part_method);
    }
    else if (dim == 3) {
      *out_pmn = PDM_generate_mesh_parallelepiped(comm,
                                                  elt_type,
                                                  1,
                                                  NULL,
                                                  -1,
                                                  -1,
                                                  -1,
                                                  2,
                                                  2,
                                                  2,
                                                  n_vtx_seg,
                                                  n_vtx_seg,
                                                  n_vtx_seg,
                                                  n_part,
                                                  part_method);
    }
    else {
      PDM_error("Invalid mesh_dimension %d (must be 2 or 3)", dim);
    }
  }
}


static void
_compute_cell_centers
(
  PDM_part_mesh_nodal_elmts_t   *pmne,
  int                            n_part,
  double                       **pvtx_coord,
  int                          **pn_pts,
  double                      ***pts_coord
)
{
  PDM_malloc(*pn_pts,    n_part, int     );
  PDM_malloc(*pts_coord, n_part, double *);

  for (int i_part = 0; i_part < n_part; i_part++) {
    int *elt_vtx_idx = NULL;
    int *elt_vtx     = NULL;
    int n_elt = PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne,
                                                               i_part,
                                                              &elt_vtx_idx,
                                                              &elt_vtx);

    (*pn_pts)[i_part] = n_elt;
    PDM_malloc((*pts_coord)[i_part], n_elt * 3, double);

    for (int i_elt = 0; i_elt < n_elt; i_elt++) {

      double *c = &(*pts_coord)[i_part][3*i_elt];

      for (int i = 0; i < 3; i++) {
        c[i] = 0;
      }

      for (int idx_vtx = elt_vtx_idx[i_elt]; idx_vtx < elt_vtx_idx[i_elt+1]; idx_vtx++) {
        int i_vtx = elt_vtx[idx_vtx] - 1;
        for (int i = 0; i < 3; i++) {
          c[i] += pvtx_coord[i_part][3*i_vtx+i];
        }
      }

      assert(elt_vtx_idx[i_elt+1] > elt_vtx_idx[i_elt]);

      double normalization = 1. / (elt_vtx_idx[i_elt+1] - elt_vtx_idx[i_elt]);
      for (int i = 0; i < 3; i++) {
        c[i] *= normalization;
      }
    }

    PDM_free(elt_vtx_idx);
    PDM_free(elt_vtx    );
  }
}




int main
(
  int   argc,
  char *argv[]
)
{
  /* Default parameters values */
  PDM_g_num_t          n_vtx_seg     = 5;
  double               length        = 1.;
  int                  n_part        = 1;
  int                  post          = 0;
  PDM_Mesh_nodal_elt_t t_elt         = PDM_MESH_NODAL_HEXA8;
  PDM_split_dual_t     part_method   = PDM_SPLIT_DUAL_WITH_HILBERT;
  char                *filename_mesh = NULL;


  /* Parse command line arguments */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
             &part_method,
             &t_elt,
             &filename_mesh);


  /* Initialize MPI */
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  /* Generate mesh */
  PDM_part_mesh_nodal_t *pmn = NULL;
  _generate_mesh(comm,
                 n_part,
                 part_method,
                 t_elt,
                 n_vtx_seg,
                 filename_mesh,
                &pmn);

  PDM_geometry_kind_t geom_kind = PDM_part_mesh_nodal_principal_geom_kind_get(pmn);

  PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);

  double **pvtx_coord = NULL;
  PDM_malloc(pvtx_coord, n_part, double *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    pvtx_coord[i_part] = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_KEEP);
  }

  /* Generate point cloud (cell centers) */
  int     *pn_elt_pts     = NULL;
  double **pelt_pts_coord = NULL;
  _compute_cell_centers(pmne,
                        n_part,
                        pvtx_coord,
                        &pn_elt_pts,
                        &pelt_pts_coord);

  /* Point location */
  int **pelt_pts_idx = NULL;
  PDM_malloc(pelt_pts_idx, n_part, int *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    pelt_pts_idx[i_part] = PDM_array_new_idx_from_const_stride_int(1, pn_elt_pts[i_part]);
  }


  double tolerance = 1e-6;

  double **distance        = NULL;
  double **projected_coord = NULL;
  int    **bary_coord_idx  = NULL;
  double **bary_coord      = NULL;
  double **uvw             = NULL;

  PDM_point_location_nodal(pmne,
                           n_part,
         (const double **) pvtx_coord,
         (const int    **) pelt_pts_idx,
         (const double **) pelt_pts_coord,
                           tolerance,
                           &distance,
                           &projected_coord,
                           &bary_coord_idx,
                           &bary_coord,
                           &uvw);


  if (post) {
    PDM_part_mesh_nodal_dump_vtk(pmn, geom_kind, "point_location_mesh");

    for (int i_part = 0; i_part < n_part; i_part++) {
      char name[99];
      sprintf(name, "point_location_proj_%d.vtk", i_rank*n_part + i_part);
      PDM_vtk_write_point_cloud(name,
                                pelt_pts_idx[i_part][pn_elt_pts[i_part]],
                                projected_coord[i_part],
                                NULL,
                                NULL);
    }
  }

  /* Finalize */
  PDM_part_mesh_nodal_free(pmn);
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(pelt_pts_idx   [i_part]);
    PDM_free(pelt_pts_coord [i_part]);
    PDM_free(distance       [i_part]);
    PDM_free(projected_coord[i_part]);
    PDM_free(bary_coord_idx [i_part]);
    PDM_free(bary_coord     [i_part]);
    PDM_free(uvw            [i_part]);
  }
  PDM_free(pn_elt_pts);
  PDM_free(pelt_pts_idx);
  PDM_free(pelt_pts_coord);
  PDM_free(distance);
  PDM_free(projected_coord);
  PDM_free(bary_coord_idx);
  PDM_free(bary_coord);
  PDM_free(uvw);
  PDM_free(pvtx_coord);

  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
