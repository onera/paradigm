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

#include "pdm_isosurface.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_generate_mesh.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


static void
_build_pmn_from_iso_result
(
  PDM_isosurface_t       *isos,
  int                     id_isosurface,
  int                     n_part,
  PDM_part_mesh_nodal_t **pmn_out,
  PDM_MPI_Comm            comm
)
{

  PDM_part_mesh_nodal_t *pmn = PDM_part_mesh_nodal_create(2, n_part, comm);
  int i_edge_section = PDM_part_mesh_nodal_section_add(pmn, PDM_MESH_NODAL_BAR2);
  // int i_face_section = PDM_part_mesh_nodal_section_add(pmn, PDM_MESH_NODAL_POLY_2D);
  
  for (int i_part=0; i_part<n_part; ++i_part) {
    // > Vertex
    double      *iso_vtx_coord = NULL;
    PDM_g_num_t *iso_vtx_gnum  = NULL;
    int iso_n_vtx = PDM_isosurface_vtx_coord_get(isos, id_isosurface, i_part, &iso_vtx_coord, PDM_OWNERSHIP_KEEP);
    PDM_isosurface_ln_to_gn_get(isos, id_isosurface, i_part, PDM_MESH_ENTITY_VTX, &iso_vtx_gnum, PDM_OWNERSHIP_KEEP);
    PDM_part_mesh_nodal_coord_set(pmn, i_part, iso_n_vtx, iso_vtx_coord, iso_vtx_gnum, PDM_OWNERSHIP_USER);

    // > Edge
    int         *iso_edge_vtx  = NULL;
    PDM_g_num_t *iso_edge_gnum = NULL;
    int iso_n_edge = PDM_isosurface_connectivity_get(isos, id_isosurface, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX, NULL, &iso_edge_vtx, PDM_OWNERSHIP_KEEP);
    PDM_isosurface_ln_to_gn_get(isos, id_isosurface, i_part, PDM_MESH_ENTITY_EDGE, &iso_edge_gnum, PDM_OWNERSHIP_KEEP);
    PDM_part_mesh_nodal_section_std_set(pmn, i_edge_section, i_part, iso_n_edge, iso_edge_vtx, iso_edge_gnum, NULL, NULL, PDM_OWNERSHIP_USER);

    // // > Face
    // int         *iso_face_vtx  = NULL;
    // PDM_g_num_t *iso_face_gnum = NULL;
    // int iso_n_face = PDM_isosurface_connectivity_get(isos, id_isosurface, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX, NULL, &iso_face_vtx, PDM_OWNERSHIP_KEEP);
    // PDM_isosurface_ln_to_gn_get(isos, id_isosurface, i_part, PDM_MESH_ENTITY_FACE, &iso_face_gnum, PDM_OWNERSHIP_KEEP);
    // PDM_part_mesh_nodal_section_std_set(pmn, i_face_section, i_part, iso_n_face, iso_face_vtx, iso_face_gnum, NULL, NULL, PDM_OWNERSHIP_USER);
  }

  // TODO: why if KEEP on pmn and USER on iso it double free ?

  // > Return
  *pmn_out = pmn;
}



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
     "  -nx      <level>  Number of vertices on the cube side (x direction).\n\n"
     "  -ny      <level>  Number of vertices on the cube side (y direction).\n\n"
     "  -l       <level>  Cube length.\n\n"
     "  -n_part  <level>  Number of partitions (if partitioned entry).\n\n"
     "  -is_dist          Is entry distributed ou partitioned.\n\n"
     "  -visu             Activate output.\n\n"
     "  -h                This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_vtx_x    Number of vertices on the cube x side
 * \param [inout]   n_vtx_y    Number of vertices on the cube y side
 * \param [inout]   n_vtx_z    Number of vertices on the cube z side
 * \param [inout]   length     Cube length
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   dist_entry Is entry distributed or partitioned (resp 1, 0)
 * \param [inout]   visu       Ensight outputs status
 *
 */
static void
_read_args(int                     argc,
           char                  **argv,
           PDM_g_num_t           *n_vtx_x,
           PDM_g_num_t           *n_vtx_y,
           PDM_g_num_t           *n_vtx_z,
           double                *length,
           int                   *n_part,
           int                   *dist_entry,
           int                   *visu)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_x = atol(argv[i]);
        *n_vtx_x = (PDM_g_num_t) _n_vtx_x;
      }
    }
    else if (strcmp(argv[i], "-ny") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_y = atol(argv[i]);
        *n_vtx_y = (PDM_g_num_t) _n_vtx_y;
      }
    }
    else if (strcmp(argv[i], "-nz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_z = atol(argv[i]);
        *n_vtx_z = (PDM_g_num_t) _n_vtx_z;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    // else if (strcmp(argv[i], "-t") == 0) {
    //   i++;
    //   if (i >= argc)
    //     _usage(EXIT_FAILURE);
    //   else {
    //     *t_elt = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    //   }
    // }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-dist_entry") == 0) {
      *dist_entry = 1;
    }
    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Init MPI
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  /*
   *  Read args
   */
  PDM_g_num_t          n_vtx_x    = 10;
  PDM_g_num_t          n_vtx_y    = 10;
  PDM_g_num_t          n_vtx_z    = 10;
  double               length     = 1.;
  int                  n_part     = 1;
  int                  order      = 1;
  PDM_Mesh_nodal_elt_t elt_type   = PDM_MESH_NODAL_TETRA4;
  int                  dist_entry = 0;
  int                  visu       = 0;
  
  _read_args(argc,
             argv,
             &n_vtx_x,
             &n_vtx_y,
             &n_vtx_z,
             &length,
             &n_part,
             &dist_entry,
             &visu);


  /*
   *  Generate mesh
   */
  PDM_dcube_nodal_t     *dcube_nodal = NULL;
  PDM_dmesh_nodal_t     *dmn         = NULL;
  PDM_part_mesh_nodal_t *pmn         = NULL;
  if (dist_entry==1) {
    dcube_nodal = PDM_dcube_nodal_gen_create(comm,
                                             n_vtx_x,
                                             n_vtx_y,
                                             n_vtx_z,
                                             length,
                                             0.,
                                             0.,
                                             0.,
                                             elt_type,
                                             order,
                                             PDM_OWNERSHIP_KEEP);
    PDM_dcube_nodal_gen_build(dcube_nodal);

    dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube_nodal);
  } else if (dist_entry==0) {
    pmn = PDM_generate_mesh_parallelepiped(comm,
                                      elt_type,
                                      order,
                                      NULL,
                                      0.,
                                      0.,
                                      0.,
                                      length,
                                      length,
                                      length,
                                      n_vtx_x,
                                      n_vtx_y,
                                      n_vtx_z,
                                      n_part,
                                      PDM_SPLIT_DUAL_WITH_PARMETIS); // TODO: Allow various partitioning ?
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_t_isosurface_2d_nodal dist_entry must be 0 or 1 (here set to %d)\n", dist_entry);
  }


  /*
   *  TODO:
   *    - extract bc
   *    - plusieurs isovalues
   *    - reequilibrate/local
   *    - test reset (et chp variable ?)
   */



  /*
   *  Creating isosurface object
   */
  PDM_isosurface_t *isos = PDM_isosurface_create(comm, 3);
  if (dist_entry==1) {
    PDM_isosurface_dmesh_nodal_set(isos, dmn);
  } else if (dist_entry==0) {
    PDM_isosurface_mesh_nodal_set(isos, pmn);
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_t_isosurface_2d_nodal dist_entry must be 0 or 1 (here set to %d)\n", dist_entry);
  }


  /*
   *  Add isosurface parameters
   */

  // > Plane isosurface
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



  /*
   *  Compute isosurface
   */
  PDM_isosurface_compute(isos, iso1);


  /*
   *  Visu isosurfaces
   */
  if (visu==1) {
    if (dist_entry==1) {
      PDM_error(__FILE__, __LINE__, 0, "PDM_t_isosurface_2d_nodal:: Not implmented\n");
    }
    else if (dist_entry==0) {
      // > iso line output
      PDM_part_mesh_nodal_t *iso_pmn1 = NULL;
      _build_pmn_from_iso_result(isos, iso1, n_part, &iso_pmn1, comm);

      PDM_part_mesh_nodal_dump_vtk(iso_pmn1,
                                   PDM_GEOMETRY_KIND_RIDGE,
                                   "pmn_iso_line_iso_mesh");

      PDM_part_mesh_nodal_free(iso_pmn1);
    }
  }



  return 0;
}
