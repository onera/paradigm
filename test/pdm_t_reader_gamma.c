
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pdm_config.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_io.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_printf.h"
#include "pdm_priv.h"
#include "pdm_reader_gamma.h"
#include "pdm_writer.h"
#include "pdm.h"


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
     "  -f <filename> Input mesh file name.\n\n"
     "  -visu         Enable export for visualization.\n\n"
     "  -h            This message.\n\n");
  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc      Number of arguments
 * \param [in]    argv      Arguments
 * \param [inout] filename  Input mesh file name
 * \param [inout] visu      Enable visualization
 *
 */
static void
_read_args
(
  int    argc,
  char **argv,
  char **filename,
  int   *visu
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



int main(int argc, char *argv[])
{
  /*
   *  Read args
   */
  char *filename = NULL;
  int   visu     = 0;
  int   n_part   = 1;

  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  _read_args(argc,
             argv,
             &filename,
             &visu);

  if (filename == NULL) {
    filename = (char *) PDM_MESH_DIR"box.mesh";
  }

  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_dmesh_nodal_t *dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                                        filename,
                                                        0,
                                                        0);

  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_domains,
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

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);

  if (visu) {

    /* Re-write in Inria format */
    PDM_part_mesh_nodal_t *pmn = NULL;
    PDM_multipart_get_part_mesh_nodal(mpart,
                                      0,
                                      &pmn,
                                      PDM_OWNERSHIP_KEEP);
    PDM_part_mesh_nodal_dump_gamma(pmn,
                                   "reader_gamma_out.mesh");

    /* Write in Ensight Gold format */
    PDM_writer_t *id_cs = PDM_writer_create("Ensight",
                                            PDM_WRITER_FMT_BIN,
                                            PDM_WRITER_TOPO_CST,
                                            PDM_WRITER_OFF,
                                            "reader_gamma",
                                            "reader_gamma",
                                            PDM_MPI_COMM_WORLD,
                                            PDM_IO_KIND_MPI_SIMPLE,
                                            1.,
                                            NULL);

    int id_geom = PDM_writer_geom_create(id_cs,
                                         "reader_gamma",
                                         n_part);

    int id_var_num_part = PDM_writer_var_create(id_cs,
                                                PDM_WRITER_ON,
                                                PDM_WRITER_VAR_SCALAR,
                                                PDM_WRITER_VAR_ELEMENTS,
                                                "num_part");

    PDM_writer_step_beg(id_cs, 0.);

    PDM_real_t **val_num_part = NULL;
    PDM_malloc(val_num_part, n_part, PDM_real_t *);

    for (int i_part = 0; i_part < n_part; i_part++) {

      int *cell_face_idx = NULL;
      int *cell_face     = NULL;
      int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &cell_face_idx,
                                                       &cell_face,
                                                       PDM_OWNERSHIP_KEEP);

      int *face_vtx_idx = NULL;
      int *face_vtx     = NULL;
      int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                       &face_vtx_idx,
                                                       &face_vtx,
                                                       PDM_OWNERSHIP_KEEP);

      double *vtx_coord = NULL;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx_coord,
                                                   PDM_OWNERSHIP_KEEP);
      PDM_g_num_t *cell_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &cell_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_g_num_t *face_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_g_num_t *vtx_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &vtx_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_writer_geom_coord_set(id_cs,
                                id_geom,
                                i_part,
                                n_vtx,
                                vtx_coord,
                                vtx_ln_to_gn,
                                PDM_OWNERSHIP_USER);

      PDM_writer_geom_cell3d_cellface_add(id_cs,
                                          id_geom,
                                          i_part,
                                          n_cell,
                                          n_face,
                                          face_vtx_idx,
                                          NULL,
                                          face_vtx,
                                          cell_face_idx,
                                          NULL,
                                          cell_face,
                                          cell_ln_to_gn);

      PDM_malloc(val_num_part[i_part], n_cell, PDM_real_t);
      for (int i = 0; i < n_cell; i++) {
        val_num_part[i_part][i] = n_part*i_rank + i_part;
      }

      PDM_writer_var_set(id_cs,
                         id_var_num_part,
                         id_geom,
                         i_part,
                         val_num_part[i_part]);
    }

    PDM_writer_geom_write(id_cs,
                          id_geom);

    PDM_writer_var_write(id_cs,
                         id_var_num_part);

    PDM_writer_var_free(id_cs,
                        id_var_num_part);

    PDM_writer_step_end(id_cs);

    for (int i = 0; i < n_part; i++) {
      PDM_free(val_num_part[i]);
    }
    PDM_free(val_num_part);

    PDM_writer_free(id_cs);
  }
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize ();

  return 0;
}
