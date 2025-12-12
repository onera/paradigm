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
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_generate_mesh.h"
#include "pdm_mem_tool.h"
#include "pdm_laplacian_smoothing.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_to_part_mesh.h"
#include "pdm_printf.h"
#include "pdm_vtk.h"
#include "pdm_part_comm_graph.h"
#include "pdm_part_comm_graph_algorithm.h"

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -mesh_path <str>     Mesh path.\n\n"
     "  -dim       <int>     Case dimension.\n\n"
     "  -n         <int>     Number of vertices per side. \n\n"
     "  -n_iter    <int>     Number of Laplacian smoothing iterations. \n\n"
     "  -damping   <double>  Damping constant of Laplacian smoothing iterations. \n\n"
     "  -tol       <double>  Tolerance for convergence of Laplacian smoothing iterations. \n\n"
     "  -visu                Enable outputs. \n\n"
     "  -h                   This message.\n\n");
  exit(exit_code);
}

static void
_read_args
(
  int                argc,
  char             **argv,
  char             **mesh_filename,
  int               *dim,
  PDM_g_num_t       *n_vtx_seg,
  PDM_split_dual_t  *part_method,
  int               *n_iter,
  double            *damping,
  double            *tol,
  int               *visu
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-in") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *mesh_filename = argv[i];
      }
    }
    else if (strcmp(argv[i], "-dim") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *dim = atoi(argv[i]);;
      }
    }
    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_vtx_seg = (PDM_g_num_t) atoi(argv[i]);;
      }
    }
    else if (strcmp(argv[i], "-part_method") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *part_method = (PDM_split_dual_t) atoi(argv[i]);;
      }
    }
    else if (strcmp(argv[i], "-n_iter") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_iter = atoi(argv[i]);;
      }
    }
    else if (strcmp(argv[i], "-damping") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *damping = atof(argv[i]);;
      }
    }
    else if (strcmp(argv[i], "-tol") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *tol = atof(argv[i]);;
      }
    }
    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}


static void
_generate_mesh
(
        PDM_MPI_Comm            comm,
        int                     n_part,
        PDM_split_dual_t        part_method,
        int                     mesh_dimension,
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
    if (mesh_dimension == 2) {
      *out_pmn = PDM_generate_mesh_rectangle(comm,
                                             PDM_MESH_NODAL_TRIA3,
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
    else if (mesh_dimension == 3) {
      *out_pmn = PDM_generate_mesh_parallelepiped(comm,
                                                  PDM_MESH_NODAL_TETRA4,
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
      PDM_error(__FILE__, __LINE__, 0, "Invalid mesh_dimension %d (must be 2 or 3)\n", mesh_dimension);
    }
  }
}


/**
 * \brief  Main
 */
int
main
(
  int   argc,
  char *argv[]
)
{
  /* Initialize MPI */
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Parse command line argmuents */
  char             *mesh_filename = NULL;
  int               dim           = 2;
  PDM_g_num_t       n_vtx_seg     = 10;
  int               n_part        = 1;
  PDM_split_dual_t  part_method   = PDM_SPLIT_DUAL_WITH_HILBERT;
  int               n_iter        = 10;
  int               visu          = 0;
  double            damping       = 1.0;
  double            tol           = 0.1;
  _read_args(argc,
             argv,
            &mesh_filename,
            &dim,
            &n_vtx_seg,
            &part_method,
            &n_iter,
            &damping,
            &tol,
            &visu);

  /* Generate mesh */
  PDM_part_mesh_nodal_t *pmn = NULL;
  _generate_mesh(comm,
                 n_part,
                 part_method,
                 dim,
                 n_vtx_seg,
                 mesh_filename,
                &pmn);

  PDM_geometry_kind_t geom_kind = PDM_part_mesh_nodal_principal_geom_kind_get(pmn);
  if (geom_kind == PDM_GEOMETRY_KIND_SURFACIC) {
    dim = 2;
  }
  else if (geom_kind == PDM_GEOMETRY_KIND_VOLUMIC) {
    dim = 3;
  }

  /* Generate edges */
  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm = PDM_part_mesh_nodal_to_part_mesh_create(pmn,
                                                                                          PDM_FALSE,
                                                                                          PDM_OWNERSHIP_USER);

  if (dim == 3) {
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_CELL_EDGE);
  }
  else {
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE);
  }
  PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX);

  PDM_part_mesh_nodal_to_part_mesh_compute(pmn_to_pm);

  PDM_part_mesh_t *pmesh = NULL;
  PDM_part_mesh_nodal_to_part_mesh_part_mesh_get(pmn_to_pm,
                                                &pmesh,
                                                 PDM_OWNERSHIP_KEEP);

  int     *pn_vtx     = NULL;
  int     *pn_edge    = NULL;
  double **pvtx_coord = NULL;
  int    **pedge_vtx  = NULL;
  PDM_malloc(pn_vtx,     n_part, int     );
  PDM_malloc(pn_edge,    n_part, int     );
  PDM_malloc(pvtx_coord, n_part, double *);
  PDM_malloc(pedge_vtx,  n_part, int    *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    pn_vtx    [i_part] = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    pvtx_coord[i_part] = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_KEEP);

    pn_edge[i_part] = PDM_part_mesh_n_entity_get(pmesh, i_part, PDM_MESH_ENTITY_EDGE);

    int *edge_vtx_idx = NULL;
    PDM_part_mesh_connectivity_get(pmesh,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   &pedge_vtx[i_part],
                                   &edge_vtx_idx,
                                   PDM_OWNERSHIP_USER);
    PDM_free(edge_vtx_idx);
  }

  PDM_part_mesh_nodal_to_part_mesh_free(pmn_to_pm);

  /* Generate field */
  double **pvtx_field = NULL;
  PDM_malloc(pvtx_field, n_part, double *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_malloc(pvtx_field[i_part], pn_vtx[i_part], double);

    for (int i_vtx = 0; i_vtx < pn_vtx[i_part]; i_vtx++) {
      double x = pvtx_coord[i_part][3*i_vtx  ];
      double y = pvtx_coord[i_part][3*i_vtx+1];
      double z = pvtx_coord[i_part][3*i_vtx+2];
      pvtx_field[i_part][i_vtx] = (x*x + y*y + z*z < 0.5) ? 1.0 : 0.0;
    }
  }

  const char    *vtx_field_name[] = {"field"};
  const double **vtx_field     [] = {(const double **) pvtx_field};

  if (visu) {
    PDM_part_mesh_nodal_dump_vtk_with_fields(pmn,
                                             geom_kind,
                                             "laplacian_init",
                                             0,
                                             NULL,
                                             NULL,
                                             1,
                                             vtx_field_name,
                                             vtx_field);
  }

  /* Generate vtx pcg */
  PDM_part_comm_graph_t *pcg_vtx = NULL;
  PDM_part_mesh_nodal_part_comm_graph_get(pmn,
                                          PDM_MESH_ENTITY_VTX,
                                          &pcg_vtx,
                                          PDM_OWNERSHIP_KEEP);

  /* Laplacian smoothing without tolerance */
  PDM_laplacian_smoothing_fields(comm,
                                 n_part,
                                 pn_vtx,
                                 pcg_vtx,
                                 NULL,
                                 NULL,
                                 pn_edge,
                                 pedge_vtx,
                                 NULL,
                                 NULL,
                                 damping,
                                 n_iter,
                                 -1.,
                                 1,
                                 pvtx_field);

  /* Generate field again */
  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i_vtx = 0; i_vtx < pn_vtx[i_part]; i_vtx++) {
      double x = pvtx_coord[i_part][3*i_vtx  ];
      double y = pvtx_coord[i_part][3*i_vtx+1];
      double z = pvtx_coord[i_part][3*i_vtx+2];
      pvtx_field[i_part][i_vtx] = (x*x + y*y + z*z < 0.5) ? 1.0 : 0.0;
    }
  }

  /* Laplacian smoothing with tolerance */
  PDM_laplacian_smoothing_fields(comm,
                                 n_part,
                                 pn_vtx,
                                 pcg_vtx,
                                 NULL,
                                 NULL,
                                 pn_edge,
                                 pedge_vtx,
                                 NULL,
                                 NULL,
                                 damping,
                                 n_iter,
                                 -1.,
                                 1,
                                 pvtx_field);

  if (visu) {
    PDM_part_mesh_nodal_dump_vtk_with_fields(pmn,
                                             geom_kind,
                                             "laplacian_final",
                                             0,
                                             NULL,
                                             NULL,
                                             1,
                                             vtx_field_name,
                                             vtx_field);
  }

  /* Generate vtx group */
  int  *pn_vtx_frozen = NULL;
  int **pvtx_frozen   = NULL;
  PDM_malloc(pn_vtx_frozen, n_part, int  );
  PDM_malloc(pvtx_frozen,   n_part, int *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    pn_vtx_frozen[i_part] = 0;
    PDM_malloc(pvtx_frozen[i_part], pn_vtx[i_part], int);
    for (int i_vtx = 0; i_vtx < pn_vtx[i_part]; i_vtx++) {
      double x = pvtx_coord[i_part][3*i_vtx];
      if (x > 0.) {
        pvtx_frozen[i_part][pn_vtx_frozen[i_part]++] = i_vtx+1;
      }
    }
    PDM_realloc(pvtx_frozen[i_part], pvtx_frozen[i_part], pn_vtx_frozen[i_part], int);
  }

  /* Generate edge pcg */
  int **pedge_vtx_idx = NULL;
  PDM_malloc(pedge_vtx_idx, n_part, int *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    pedge_vtx_idx[i_part] = PDM_array_new_idx_from_const_stride_int(2, pn_edge[i_part]);
  }
  PDM_part_comm_graph_t *pcg_edge = NULL;
  PDM_part_comm_graph_entity1_to_part_comm_graph_entity2(pcg_vtx,
                                                         pn_vtx,
                                                         pn_edge,
                                                         pedge_vtx_idx,
                                                         pedge_vtx,
                                                         &pcg_edge);
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(pedge_vtx_idx[i_part]);
  }
  PDM_free(pedge_vtx_idx);

  /* Generate edge weights */
  double **pedge_weight = NULL;
  PDM_laplacian_smoothing_idw_weights_compute(n_part,
                                              pvtx_coord,
                                              pn_edge,
                                              pedge_vtx,
                                              2,
                                              &pedge_weight);

  /* Generate strided fields */
  double **pvtx_coord_tmp = NULL;
  PDM_malloc(pvtx_coord_tmp, n_part, double *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_malloc(pvtx_coord_tmp[i_part], 3 * pn_vtx[i_part], double);
    memcpy(pvtx_coord_tmp[i_part], pvtx_coord[i_part], sizeof(double) * 3 * pn_vtx[i_part]);
  }

  /* Laplacian smoothing with tolerance, groups, pcg_edge, weights and strided fields */
  PDM_laplacian_smoothing_fields(comm,
                                 n_part,
                                 pn_vtx,
                                 pcg_vtx,
                                 pn_vtx_frozen,
                                 pvtx_frozen,
                                 pn_edge,
                                 pedge_vtx,
                                 pedge_weight,
                                 pcg_edge,
                                 damping,
                                 n_iter,
                                 -1.,
                                 3,
                                 pvtx_coord_tmp);

  /* Finalize */
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(pedge_vtx   [i_part]);
    PDM_free(pvtx_field  [i_part]);
    PDM_free(pvtx_frozen [i_part]);
    PDM_free(pedge_weight[i_part]);
    PDM_free(pvtx_coord_tmp[i_part]);
  }
  PDM_free(pn_vtx    );
  PDM_free(pn_edge   );
  PDM_free(pvtx_coord);
  PDM_free(pedge_vtx );
  PDM_free(pvtx_field);
  PDM_free(pvtx_frozen   );
  PDM_free(pn_vtx_frozen );
  PDM_free(pedge_weight  );
  PDM_free(pvtx_coord_tmp);
  PDM_part_mesh_nodal_free(pmn);
  PDM_part_comm_graph_free(pcg_edge);

  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
