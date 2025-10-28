#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_generate_mesh.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_part_comm_graph.h"
#include "pdm_printf.h"
#include "pdm_vtk.h"
#include "pdm.h"


static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -mesh_path  Mesh path.\n\n"
     "  -dim        Mesh dimension.\n\n"
     "  -visu       Enable outputs.\n\n"
     "  -h          This message.\n\n");
  exit(exit_code);
}


static void
_read_args
(
  int     argc,
  char  **argv,
  char  **mesh_path,
  int    *dim,
  int    *visu
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-mesh_path") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *mesh_path = argv[i];
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
  char             *mesh_path   = NULL;
  int               dim         = 3;
  PDM_g_num_t       n_vtx_seg   = 10;
  int               n_part      = 1;
  PDM_split_dual_t  part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  int               visu        = 0;
  _read_args(argc,
             argv,
            &mesh_path,
            &dim,
            &visu);


  /* Generate mesh */
  PDM_part_mesh_nodal_t *pmn = NULL;
  _generate_mesh(comm,
                 n_part,
                 part_method,
                 dim,
                 n_vtx_seg,
                 mesh_path,
                &pmn);

  /* Select elements in group 1 of dimension dim-1 */
  PDM_geometry_kind_t geom_kind = PDM_part_mesh_nodal_principal_geom_kind_get(pmn) + 1;
  int i_group = PDM_part_mesh_nodal_n_group_get(pmn, geom_kind) - 1;

  int  *n_selected_elt = NULL;
  int **selected_elt   = NULL;
  PDM_malloc(n_selected_elt, n_part, int  );
  PDM_malloc(selected_elt,   n_part, int *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_g_num_t *group_ln_to_gn = NULL;
    PDM_part_mesh_nodal_group_get(pmn,
                                  geom_kind,
                                  i_part,
                                  i_group,
                                  &n_selected_elt[i_part],
                                  &selected_elt  [i_part],
                                  &group_ln_to_gn,
                                  PDM_OWNERSHIP_BAD_VALUE);
  }

  /* Get element->vertex connectivity */
  int **elt_vtx_idx = NULL;
  int **elt_vtx     = NULL;
  PDM_malloc(elt_vtx_idx, n_part, int *);
  PDM_malloc(elt_vtx,     n_part, int *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_part_mesh_nodal_cell_vtx_connect_get(pmn,
                                             geom_kind,
                                             i_part,
                                            &elt_vtx_idx[i_part],
                                            &elt_vtx    [i_part]);
  }

  /* Get vertex Part Comm Graph */
  PDM_part_comm_graph_t *pcg_vtx = NULL;
  PDM_part_mesh_nodal_part_comm_graph_get(pmn,
                                          PDM_MESH_ENTITY_VTX,
                                          &pcg_vtx,
                                          PDM_OWNERSHIP_BAD_VALUE);

  /* Deduce selected vertices */
  int  *n_selected_vtx = NULL;
  int **selected_vtx   = NULL;
  PDM_part_comm_graph_selected_entity1_to_selected_entity2(n_selected_elt,
                                                           selected_elt,
                                                           elt_vtx_idx,
                                                           elt_vtx,
                                                           pcg_vtx,
                                                          &n_selected_vtx,
                                                          &selected_vtx);

  /* Visualize */
  if (visu) {
    PDM_part_mesh_nodal_dump_vtk(pmn, geom_kind-1, "pcg_group_vtx_mesh");
    PDM_part_mesh_nodal_dump_vtk(pmn, geom_kind,   "pcg_group_vtx_elt");

    for (int i_part = 0; i_part < n_part; i_part++) {

      double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);

      double *coord = NULL;
      PDM_malloc(coord, n_selected_vtx[i_part] * 3, double);
      for (int i = 0; i < n_selected_vtx[i_part]; i++) {
        int i_vtx = selected_vtx[i_part][i] - 1;
        memcpy(&coord[3*i], &vtx_coord[3*i_vtx], sizeof(double) * 3);
      }

      char name[99];
      sprintf(name, "pcg_group_vtx_%d.vtk", i_rank*n_part+i_part);
      PDM_vtk_write_point_cloud(name,
                                n_selected_vtx[i_part],
                                coord,
                                NULL,
                                selected_vtx[i_part]);
      PDM_free(coord);
    }
  }

  /* Finalize */
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(elt_vtx_idx [i_part]);
    PDM_free(elt_vtx     [i_part]);
    PDM_free(selected_vtx[i_part]);
  }
  PDM_free(elt_vtx_idx   );
  PDM_free(elt_vtx       );
  PDM_free(n_selected_vtx);
  PDM_free(  selected_vtx);
  PDM_free(n_selected_elt);
  PDM_free(  selected_elt);
  PDM_part_mesh_nodal_free(pmn);

  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}