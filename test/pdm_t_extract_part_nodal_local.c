
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pdm.h"
#include "pdm_array.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_extract_part.h"
#include "pdm_mem_tool.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_reader_gamma.h"

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
_usage
(
  int exit_code
)
{
  printf
    ("\n"
     "  -h                      This message.\n\n"
     "  -n_part     <n>         Number of partitions.\n\n"
     "  -elt_type   <t>         Element type.\n\n"
     "  -in         <filename>  Mesh file name.\n\n"
     "  -visu                   Enable exports for visualization.\n\n"
     "  -n          <n>         Number of vertices per side (only for automatically generated mesh).\n\n"
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
  PDM_Mesh_nodal_elt_t  *elt_type,
  char                 **filename,
  PDM_split_dual_t      *part_method,
  int                   *visu,
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

    else if (strcmp(argv[i], "-elt_type") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-in") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *filename = argv[i];
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
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

    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
    }

    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    }

    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}


static
PDM_part_mesh_nodal_t *
_gen_mesh
(
 PDM_MPI_Comm          comm,
 const char           *mesh_file,
 int                   n_part,
 PDM_g_num_t           n_vtx_seg,
 PDM_Mesh_nodal_elt_t  elt_type,
 PDM_split_dual_t      part_method
)
{
  PDM_dmesh_nodal_t *dmn   = NULL;
  PDM_dmesh_t       *dmesh = NULL;

  int n_domain = 1;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  if (mesh_file == NULL) {
    if (elt_type == PDM_MESH_NODAL_POLY_3D) {
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
                       0,
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
      dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
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

      PDM_multipart_dmesh_set(mpart, 0, dmesh);
    }
    else {
      PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                            n_vtx_seg,
                                                            n_vtx_seg,
                                                            n_vtx_seg,
                                                            1.,
                                                            0.,
                                                            0.,
                                                            0.,
                                                            elt_type,
                                                            1,
                                                            PDM_OWNERSHIP_USER);
      PDM_dcube_nodal_gen_build(dcube);

      dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
      PDM_dcube_nodal_gen_free(dcube);

      PDM_dmesh_nodal_generate_distribution(dmn);

      PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
    }
  }
  else {
    dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                       mesh_file,
                                       0,
                                       0);
    PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  }


  PDM_multipart_compute(mpart);


  /* Let's go */
  PDM_part_mesh_nodal_t *pmn = NULL;

  if (elt_type == PDM_MESH_NODAL_POLY_3D) {
    PDM_dmesh_free(dmesh);

    pmn = PDM_part_mesh_nodal_create(3, n_part, comm);

    for (int i_part = 0; i_part < n_part; i_part++) {

      int *cell_face_idx = NULL;
      int *cell_face     = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                          &cell_face_idx,
                                          &cell_face,
                                          PDM_OWNERSHIP_KEEP);

      int *face_vtx_idx = NULL;
      int *face_vtx     = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          &face_vtx_idx,
                                          &face_vtx,
                                          PDM_OWNERSHIP_KEEP);

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

      PDM_g_num_t *vtx_ln_to_gn = NULL;
      int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                  0,
                                                  i_part,
                                                  PDM_MESH_ENTITY_VTX,
                                                  &vtx_ln_to_gn,
                                                  PDM_OWNERSHIP_USER);

      double *vtx_coord = NULL;
      PDM_multipart_part_vtx_coord_get(mpart,
                                       0,
                                       i_part,
                                       &vtx_coord,
                                       PDM_OWNERSHIP_USER);

      PDM_part_mesh_nodal_coord_set(pmn,
                                    i_part,
                                    n_vtx,
                                    vtx_coord,
                                    PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_nodal_vtx_gnum_set(pmn,
                                       i_part,
                                       vtx_ln_to_gn,
                                       PDM_OWNERSHIP_KEEP);

      PDM_part_mesh_nodal_cell3d_cellface_add(pmn,
                                              i_part,
                                              n_cell,
                                              n_face,
                                              face_vtx_idx,
                                              face_vtx,
                                              face_ln_to_gn,
                                              cell_face_idx,
                                              cell_face,
                                              cell_ln_to_gn,
                                              PDM_OWNERSHIP_KEEP);
    }
  }
  else {
    PDM_multipart_get_part_mesh_nodal(mpart,
                                      0,
                                      &pmn,
                                      PDM_OWNERSHIP_USER);
  }

  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);


  PDM_geometry_kind_t geom_kind = PDM_part_mesh_nodal_principal_geom_kind_get(pmn);

  if (geom_kind == PDM_GEOMETRY_KIND_VOLUMIC) {
    PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn,
                                                                                      geom_kind);
    int n_group = PDM_part_mesh_nodal_elmts_n_group_get(pmne);

    if (n_group == 0) {
      // Add volume groups
      n_group = 3;

      PDM_part_mesh_nodal_elmts_n_group_set(pmne, n_group);

      int         **group_elt          = NULL;
      PDM_g_num_t **group_elt_ln_to_gn = NULL;
      PDM_malloc(group_elt,          n_group, int         *);
      PDM_malloc(group_elt_ln_to_gn, n_group, PDM_g_num_t *);

      for (int i_part = 0; i_part < n_part; i_part++) {
        int n_elt = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);

        PDM_g_num_t *g_num = PDM_part_mesh_nodal_elmts_g_num_get_from_part(pmne,
                                                                           i_part,
                                                                           PDM_OWNERSHIP_KEEP);
        int *group_elt_n = PDM_array_zeros_int(n_group);
        for (int i_elt = 0; i_elt < n_elt; i_elt++) {
          int i_group = g_num[i_elt] % n_group;
          group_elt_n[i_group]++;
        }

        for (int i_group = 0; i_group < n_group; i_group++) {
          PDM_malloc(group_elt         [i_group], group_elt_n[i_group], int        );
          PDM_malloc(group_elt_ln_to_gn[i_group], group_elt_n[i_group], PDM_g_num_t);
          group_elt_n[i_group] = 0;
        }

        for (int i_elt = 0; i_elt < n_elt; i_elt++) {
          int i_group = g_num[i_elt] % n_group;
          group_elt         [i_group][group_elt_n[i_group]] = i_elt+1;
          group_elt_ln_to_gn[i_group][group_elt_n[i_group]] = 1 + (g_num[i_elt] - 1) / n_group;
          group_elt_n[i_group]++;
        }

        for (int i_group = 0; i_group < n_group; i_group++) {
          PDM_part_mesh_nodal_elmts_group_set(pmne,
                                              i_part,
                                              i_group,
                                              group_elt_n       [i_group],
                                              group_elt         [i_group],
                                              group_elt_ln_to_gn[i_group],
                                              PDM_OWNERSHIP_KEEP);
        }
        PDM_free(group_elt_n);
      }
      PDM_free(group_elt         );
      PDM_free(group_elt_ln_to_gn);
    }
  }

  return pmn;
}


static void
_get_extracted
(
 PDM_MPI_Comm                   comm,
 PDM_part_mesh_nodal_elmts_t   *pmne,
 int                            n_part,
 int                          **out_n_extract,
 int                         ***out_extract_lnum
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  int  *n_extract    = PDM_array_zeros_int(n_part);
  int **extract_lnum = NULL;
  PDM_malloc(extract_lnum, n_part, int *);

  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i_section = 0; i_section < n_section; i_section++) {
      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, sections_id[i_section], i_part);

      PDM_g_num_t *ln_to_gn = PDM_part_mesh_nodal_elmts_g_num_get(pmne,
                                                                  sections_id[i_section],
                                                                  i_part,
                                                                  PDM_OWNERSHIP_KEEP);
      for (int i_elt = 0; i_elt < n_elt; i_elt++) {
        if (ln_to_gn[i_elt]%5 == 0) {
          n_extract[i_part]++;
        }
      }
    }


    PDM_malloc(extract_lnum[i_part], n_extract[i_part], int);

    n_extract[i_part] = 0;


    int i_cell = -1;
    for (int i_section = 0; i_section < n_section; i_section++) {

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, sections_id[i_section], i_part);


      PDM_g_num_t *ln_to_gn = PDM_part_mesh_nodal_elmts_g_num_get(pmne,
                                                                  sections_id[i_section],
                                                                  i_part,
                                                                  PDM_OWNERSHIP_KEEP);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 sections_id[i_section],
                                                                 i_part,
                                                                 PDM_OWNERSHIP_KEEP);

      for (int i_elt = 0; i_elt < n_elt; i_elt++) {
        if (ln_to_gn[i_elt]%5 == 0) {
          if (parent_num == NULL) {
            i_cell++;
          }
          else {
            i_cell = parent_num[i_elt];
          }
          extract_lnum[i_part][n_extract[i_part]++] = i_cell + 1;
        }
      }
    }

    // PDM_log_trace_array_int(extract_lnum[i_part], n_extract[i_part], "extract_lnum : ");
  }

  *out_n_extract    = n_extract;
  *out_extract_lnum = extract_lnum;
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
  int                   n_part      = 1;
  PDM_Mesh_nodal_elt_t  elt_type    = PDM_MESH_NODAL_PRISM6;
  char                 *mesh_file   = NULL;
  PDM_split_dual_t      part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  int                   visu        = 0;
  PDM_g_num_t           n_vtx_seg   = 10;

  PDM_bool_t compute_child_gnum = PDM_TRUE;

  _read_args(argc,
             argv,
             &n_part,
             &elt_type,
             &mesh_file,
             &part_method,
             &visu,
             &n_vtx_seg);

  /* Generate mesh */
  PDM_part_mesh_nodal_t *pmn = _gen_mesh(comm,
                                         mesh_file,
                                         n_part,
                                         n_vtx_seg,
                                         elt_type,
                                         part_method);

  int mesh_dimension = PDM_part_mesh_nodal_mesh_dimension_get(pmn);

  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
  if (mesh_dimension == 2) {
    geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
  }
  else if (mesh_dimension == 1) {
    geom_kind = PDM_GEOMETRY_KIND_RIDGE;
  }

  if (visu) {
    for (PDM_geometry_kind_t geom_kind_child = geom_kind; geom_kind_child < PDM_GEOMETRY_KIND_CORNER; geom_kind_child++) {
      char name[999];
      sprintf(name, "init_pmn_%d", geom_kind_child);
      PDM_part_mesh_nodal_dump_vtk(pmn,
                                   geom_kind_child,
                                   name);
    }
  }



  /* Selection */
  PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn,
                                                                                    geom_kind);
  int  *n_extract    = NULL;
  int **extract_lnum = NULL;
  _get_extracted(comm,
                 pmne,
                 n_part,
                 &n_extract,
                 &extract_lnum);


  /* Extract part */
  PDM_extract_part_t *extrp = PDM_extract_part_create(mesh_dimension,
                                                      n_part,
                                                      n_part,
                                                      PDM_EXTRACT_PART_KIND_LOCAL,
                                                      part_method,
                                                      compute_child_gnum,
                                                      PDM_OWNERSHIP_KEEP,
                                                      comm);

  PDM_extract_part_part_nodal_set(extrp, pmn);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_extract_part_selected_lnum_set(extrp,
                                       i_part,
                                       n_extract   [i_part],
                                       extract_lnum[i_part],
                                       PDM_OWNERSHIP_USER);
  }

  PDM_extract_part_compute(extrp);


  if (visu) {
    PDM_part_mesh_nodal_t *extract_pmn = NULL;
    PDM_extract_part_part_mesh_nodal_get(extrp,
                                         &extract_pmn,
                                         PDM_OWNERSHIP_KEEP);

    for (PDM_geometry_kind_t geom_kind_child = geom_kind; geom_kind_child < PDM_GEOMETRY_KIND_CORNER; geom_kind_child++) {
      char name[999];
      sprintf(name, "extract_pmn_%d", geom_kind_child);
      PDM_part_mesh_nodal_dump_vtk(extract_pmn,
                                   geom_kind_child,
                                   name);
    }
  }


  /* Free memory */
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(extract_lnum[i_part]);
  }
  PDM_free(extract_lnum);
  PDM_free(n_extract   );

  PDM_extract_part_free(extrp);
  PDM_part_mesh_nodal_free(pmn);

  if (i_rank == 0) {
    printf("End :D\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
