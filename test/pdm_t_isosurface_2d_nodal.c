#include <stdlib.h>
#include <string.h>

#include "pdm.h"
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
  double               length     = 1.;
  int                  n_part     = 1;
  int                  order      = 1;
  PDM_Mesh_nodal_elt_t elt_type   = PDM_MESH_NODAL_TRIA3;
  int                  dist_entry = 0;
  int                  visu       = 0;
  
  _read_args(argc,
             argv,
             &n_vtx_x,
             &n_vtx_y,
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
  log_trace("\n\n");
  log_trace("===============\n");
  log_trace("> Generate mesh\n");
  if (dist_entry==1) {
    dcube_nodal = PDM_dcube_nodal_gen_create(comm,
                                             n_vtx_x,
                                             n_vtx_y,
                                             0,
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
    pmn  = PDM_generate_mesh_rectangle(comm,
                                       elt_type,
                                       order,
                                       NULL,
                                       0.,
                                       0.,
                                       0.,
                                       length,
                                       length,
                                       n_vtx_x,
                                       n_vtx_y,
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
  log_trace("\n\n");
  log_trace("============================\n");
  log_trace("> Creating isosurface object\n");
  PDM_isosurface_t *isos = PDM_isosurface_create(comm,
                                                 2/*,
                                                 PDM_MESH_NODAL_BAR2*/);
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
  double plane_isovalues[3] = {-0.25,0.,0.25};
  int iso1 = PDM_isosurface_add(isos, 
                                PDM_ISO_SURFACE_KIND_PLANE,
                                3,
                                plane_isovalues);
  PDM_isosurface_equation_set(isos,
                              iso1,
                              plane_equation,
                              0);

  // > Sphere isosurface
  double sphere_equation [4] = {0.5,0.5,0.5,0.5};
  double sphere_isovalues[2] = {0.,0.1};
  int iso2 = PDM_isosurface_add(isos, 
                                PDM_ISO_SURFACE_KIND_SPHERE,
                                2,
                                sphere_isovalues);
  PDM_isosurface_equation_set(isos,
                              iso2,
                              sphere_equation,
                              0);

  // > User field isosurface
  double  *dfield = NULL;
  double **field  = NULL;

  double field_isovalues[1] = {0.};
  int iso3 = PDM_isosurface_add(isos, 
                                PDM_ISO_SURFACE_KIND_FIELD,
                                1,
                                field_isovalues);
  if (dist_entry==1) {
    PDM_isosurface_dfield_set(isos,
                              iso3,
                              dfield);
  } else if (dist_entry==0) {
    for (int i_part=0; i_part<n_part; ++i_part) {
      PDM_isosurface_field_set(isos,
                               iso3,
                               i_part,
                               field[i_part]);
    }
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_t_isosurface_2d_nodal dist_entry must be 0 or 1 (here set to %d)\n", dist_entry);
  }


  PDM_isosurface_compute(isos, iso1);

  /*
   *  Free objects
   */
  PDM_isosurface_free(isos);
  if (dist_entry==1) {
    PDM_dcube_nodal_gen_free(dcube_nodal);
  } else if (dist_entry==0) {
    PDM_part_mesh_nodal_free(pmn);
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_t_isosurface_2d_nodal dist_entry must be 0 or 1 (here set to %d)\n", dist_entry);
  }


  PDM_MPI_Finalize();

  return 0;
}
