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



static const int it_max = 40;
static void
_eval_mandelbrot
(
 const double  x,
 const double  y,
 const double  z,
       double *f
)
{
  PDM_UNUSED(z);

  double _x = x - 0.5;
  double _y = y;

  double xk = 0;
  double yk = 0;

  double zxk = 1.;
  double zyk = 0.;
  double rk;

  int it;
  for (it = 0; it < it_max; it++) {
    double xk_new = xk*xk - yk*yk + _x;
    double yk_new = 2*xk*yk       + _y;
    xk = xk_new;
    yk = yk_new;
    rk = sqrt(xk*xk + yk*yk);

    double zxk_new = 2*(xk*zxk - yk*zyk) + 1;
    double zyk_new = 2*(zxk*yk + xk*zyk) + 1;

    zxk = zxk_new;
    zyk = zyk_new;

    if (rk > 2.) {
      break;
    }
  }

  double mag_dz = sqrt(zxk*zxk + zyk*zyk);
  *f = rk * log(PDM_MAX(1e-9, rk)) / PDM_MAX(1e-9, mag_dz);
  // if (it < it_max-1) {
  //   *f = 1.;
  // } else {
  //   *f = -1.;
  // }
}



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

  PDM_part_mesh_nodal_t *pmn = PDM_part_mesh_nodal_create(1, n_part, comm);
  int i_edge_section = PDM_part_mesh_nodal_section_add(pmn, PDM_MESH_NODAL_BAR2);
  
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
     "  -nx       <n>  Number of vertices on the cube side (x direction).\n\n"
     "  -ny       <n>  Number of vertices on the cube side (y direction).\n\n"
     "  -l        <n>  Cube length.\n\n"
     "  -elt_type <t>  Surface element type (only for automatically generated mesh).\n\n"
     "  -n_part   <n>  Number of partitions (if partitioned entry).\n\n"
     "  -is_dist       Is entry distributed ou partitioned.\n\n"
     "  -local         Deactivate isosurface redistribution.\n\n"
     "  -visu          Activate output.\n\n"
     "  -h             This message.\n\n");

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
 * \param [inout]   length     Square length
 * \param [inout]   elt_type   Element type
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   dist_entry Is entry distributed or partitioned (resp 1, 0)
 * \param [inout]   local      Deactivate isosurface redistribution
 * \param [inout]   visu       Ensight outputs status
 *
 */
static void
_read_args(int                     argc,
           char                  **argv,
           PDM_g_num_t           *n_vtx_x,
           PDM_g_num_t           *n_vtx_y,
           double                *length,
           PDM_Mesh_nodal_elt_t  *elt_type,
           int                   *n_part,
           int                   *dist_entry,
           int                   *local,
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
    else if (strcmp(argv[i], "-elt_type") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
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
    else if (strcmp(argv[i], "-local") == 0) {
      *local = 1;
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
  PDM_g_num_t          n_vtx_x    = 5;
  PDM_g_num_t          n_vtx_y    = 5;
  double               length     = 1.;
  int                  n_part     = 1;
  int                  order      = 1;
  PDM_Mesh_nodal_elt_t elt_type   = PDM_MESH_NODAL_TRIA3;
  int                  dist_entry = 0;
  int                  local      = 0;
  int                  visu       = 0;
  
  _read_args(argc,
             argv,
             &n_vtx_x,
             &n_vtx_y,
             &length,
             &elt_type,
             &n_part,
             &dist_entry,
             &local,
             &visu);

  if (n_part <= 0) {
    dist_entry = 1;
  }


  /*
   *  Generate mesh
   */
  PDM_dcube_nodal_t     *dcube_nodal = NULL;
  PDM_dmesh_nodal_t     *dmn         = NULL;
  PDM_part_mesh_nodal_t *pmn         = NULL;
  if (dist_entry==1) {
    // Block-distributed
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
  }
  else if (dist_entry==0) {
    // Partitioned
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
  PDM_isosurface_t *isos = PDM_isosurface_create(comm,
                                                 2/*,
                                                 PDM_MESH_NODAL_BAR2*/);
  if (dist_entry==1) {
    PDM_isosurface_dmesh_nodal_set(isos, dmn);
  } else if (dist_entry==0) {
    PDM_isosurface_mesh_nodal_set(isos, pmn);
    if (local==0) {
      PDM_isosurface_redistribution_set(isos, PDM_EXTRACT_PART_KIND_REEQUILIBRATE, PDM_SPLIT_DUAL_WITH_HILBERT); // TODO: Test various partitioning ?
    }
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
    int     dn_vtx     = PDM_DMesh_nodal_n_vtx_get(dmn);
    double *dvtx_coord = PDM_DMesh_nodal_vtx_get  (dmn);
    PDM_malloc(dfield, dn_vtx, double);
    for (int i_vtx=0; i_vtx<dn_vtx; ++i_vtx) {
      _eval_mandelbrot(dvtx_coord[3*i_vtx  ],
                       dvtx_coord[3*i_vtx+1],
                       dvtx_coord[3*i_vtx+2],
                       &dfield[i_vtx]);
    }
    PDM_isosurface_dfield_set(isos,
                              iso3,
                              dfield);
  } else if (dist_entry==0) {
    PDM_malloc(field, n_part, double *);

    for (int i_part=0; i_part<n_part; ++i_part) {
      int     n_vtx     = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);
      PDM_malloc(field[i_part], n_vtx, double);
      
      for (int i_vtx=0; i_vtx<n_vtx; ++i_vtx) {
        _eval_mandelbrot( vtx_coord[3*i_vtx  ],
                          vtx_coord[3*i_vtx+1],
                          vtx_coord[3*i_vtx+2],
                         &field[i_part][i_vtx]);
      }
      PDM_isosurface_field_set(isos,
                               iso3,
                               i_part,
                               field[i_part]);
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_t_isosurface_2d_nodal dist_entry must be 0 or 1 (here set to %d)\n", dist_entry);
  }


  /*
   *  Compute isosurfaces
   */
  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_reset  (isos, iso1);
  PDM_isosurface_reset  (isos, iso2);
  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_compute(isos, iso2);
  PDM_isosurface_compute(isos, iso3);

  /*
   *  Visu isosurfaces
   */
  if (visu==1) {
    if (dist_entry==1) {
      PDM_free(dfield);
      // PDM_error(__FILE__, __LINE__, 0, "PDM_t_isosurface_2d_nodal:: Not implmented\n");
    }
    else if (dist_entry==0) {
      // > iso line output
      PDM_part_mesh_nodal_t *iso_pmn1 = NULL;
      _build_pmn_from_iso_result(isos, iso1, n_part, &iso_pmn1, comm);

      PDM_part_mesh_nodal_dump_vtk(iso_pmn1,
                                   PDM_GEOMETRY_KIND_RIDGE,
                                   "pmn_iso_line_iso_mesh");

      PDM_part_mesh_nodal_free(iso_pmn1);


      // > iso circle output
      PDM_part_mesh_nodal_t *iso_pmn2 = NULL;
      _build_pmn_from_iso_result(isos, iso2, n_part, &iso_pmn2, comm);

      PDM_part_mesh_nodal_dump_vtk(iso_pmn2,
                                   PDM_GEOMETRY_KIND_RIDGE,
                                   "pmn_iso_circle_mesh");

      PDM_part_mesh_nodal_free(iso_pmn2);


      // > iso field output
      PDM_part_mesh_nodal_t *iso_pmn3 = NULL;
      _build_pmn_from_iso_result(isos, iso3, n_part, &iso_pmn3, comm);

      PDM_part_mesh_nodal_dump_vtk(iso_pmn3,
                                   PDM_GEOMETRY_KIND_RIDGE,
                                   "pmn_iso_field_mesh");

      PDM_part_mesh_nodal_free(iso_pmn3);
    }
  }


  /*
   * TODO: 
   *   - if default test conf, put assert ?
   */


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

  if (field!=NULL) {
    for (int i_part=0; i_part<n_part; ++i_part) {
      if (field[i_part]!=NULL) {
        PDM_free(field[i_part]);
      }
    }
    PDM_free(field);
  }


  PDM_MPI_Finalize();

  return 0;
}
