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

#include "pdm_isosurface.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_generate_mesh.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


typedef double (*signed_dist_func_t)
(
  double  x,
  double  y,
  double  z
);

// Quaternion multiplication
static void
_quat_mult
(
  const double a[4],
  const double b[4],
  double       ab[4]
)
{
  ab[0] = a[0]*b[0] - PDM_DOT_PRODUCT(&a[1], &b[1]);

  PDM_CROSS_PRODUCT(&ab[1], &a[1], &b[1]);
  for (int i = 0; i < 3; i++) {
    ab[i+1] += a[0]*b[i+1] + b[0]*a[i+1];
  }
}

static int
_julia_iterations
(
  double q[4],
  double dq[4],
  const double c[4],
  const int max_iter
)
{
  for (int iter = 0; iter < max_iter; iter++) {
    double tmp[4];

    _quat_mult(q, dq, tmp);
    for (int i = 0; i < 4; i++) {
      dq[i] = 2*tmp[i];
    }

    double mag2 = 0;
    _quat_mult(q, q, tmp);
    for (int i = 0; i < 4; i++) {
      q[i] = c[i] + tmp[i];

      mag2 += q[i]*q[i];
    }

    if (mag2 > 1.e4) {//4) {
      return iter+1;
    }
  }

  return 0;
}


static double
_eval_julia4d
(
  const double x,
  const double y,
  const double z,
  const double c[4],
  const int return_it
)
{
  const int max_iter = 10;//20;
  double  q[4] = {x, y, z, 0};
  double dq[4] = {1, 0, 0, 0};

  int stat = _julia_iterations(q, dq, c, max_iter);
  PDM_UNUSED(stat);

  double mag_q  = 0;
  double mag_dq = 0;
  for (int i = 0; i < 4; i++) {
    mag_q  +=  q[i] *  q[i];
    mag_dq += dq[i] * dq[i];
  }

  mag_q  = sqrt(mag_q);
  mag_dq = sqrt(mag_dq);
  double dist = 0.5 * mag_q * log(mag_q) / mag_dq;

  if (dist < 0) {
    dist = atan(0.01*dist);
  }

  if (return_it==1) {
    return (double) stat;
  }
  else {
    return dist;
  }
}


static void
_build_pmn_from_iso_result
(
  PDM_isosurface_t       *isos,
  int                     id_iso,
  int                     n_part,
  PDM_part_mesh_nodal_t **pmn_out,
  PDM_MPI_Comm            comm
)
{

  PDM_part_mesh_nodal_t *pmn = PDM_part_mesh_nodal_create(2, n_part, comm);
  int i_edge_section = PDM_part_mesh_nodal_section_add(pmn, PDM_MESH_NODAL_BAR2);
  int i_face_section = PDM_part_mesh_nodal_section_add(pmn, PDM_MESH_NODAL_POLY_2D);

  for (int i_part=0; i_part<n_part; ++i_part) {
    // > Vertex
    double      *iso_vtx_coord = NULL;
    PDM_g_num_t *iso_vtx_gnum  = NULL;
    int iso_n_vtx = PDM_isosurface_vtx_coord_get(isos, id_iso, i_part, &iso_vtx_coord, PDM_OWNERSHIP_KEEP);
    PDM_isosurface_ln_to_gn_get(isos, id_iso, i_part, PDM_MESH_ENTITY_VTX, &iso_vtx_gnum, PDM_OWNERSHIP_KEEP);
    PDM_part_mesh_nodal_coord_set(pmn, i_part, iso_n_vtx, iso_vtx_coord, iso_vtx_gnum, PDM_OWNERSHIP_USER);

    // > Edge
    int         *iso_edge_vtx  = NULL;
    PDM_g_num_t *iso_edge_gnum = NULL;
    int iso_n_edge = PDM_isosurface_connectivity_get(isos, id_iso, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX, NULL, &iso_edge_vtx, PDM_OWNERSHIP_KEEP);
    PDM_isosurface_ln_to_gn_get(isos, id_iso, i_part, PDM_MESH_ENTITY_EDGE, &iso_edge_gnum, PDM_OWNERSHIP_KEEP);
    PDM_part_mesh_nodal_section_std_set(pmn, i_edge_section, i_part, iso_n_edge, iso_edge_vtx, iso_edge_gnum, NULL, NULL, PDM_OWNERSHIP_USER);

    // > Face
    int         *iso_face_vtx_idx  = NULL;
    int         *iso_face_vtx      = NULL;
    PDM_g_num_t *iso_face_gnum     = NULL;
    int iso_n_face = PDM_isosurface_connectivity_get(isos, id_iso, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX, &iso_face_vtx_idx, &iso_face_vtx, PDM_OWNERSHIP_KEEP);
    PDM_isosurface_ln_to_gn_get(isos, id_iso, i_part, PDM_MESH_ENTITY_FACE, &iso_face_gnum, PDM_OWNERSHIP_KEEP);
    PDM_part_mesh_nodal_section_poly2d_set(pmn, i_face_section, i_part, iso_n_face, iso_face_vtx_idx, iso_face_vtx, iso_face_gnum, NULL, PDM_OWNERSHIP_USER);
  }

  // TODO: why if KEEP on pmn and USER on iso it double free ?

  // > Return
  *pmn_out = pmn;
}


static void
_output_iso_result
(
  PDM_isosurface_t       *isos,
  int                     id_iso,
  int                     n_part,
  double               ***iso_vtx_field,
  PDM_MPI_Comm            comm
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  
  for (int i_part=0; i_part<n_part; ++i_part) {

    /**
     * Get isosurface data
     */
    // > Vertex
    double      *iso_vtx_coord = NULL;
    PDM_g_num_t *iso_vtx_gnum  = NULL;
    int iso_n_vtx = PDM_isosurface_vtx_coord_get(isos, id_iso, i_part, &iso_vtx_coord, PDM_OWNERSHIP_KEEP);
    PDM_isosurface_ln_to_gn_get(isos, id_iso, i_part, PDM_MESH_ENTITY_VTX, &iso_vtx_gnum, PDM_OWNERSHIP_KEEP);

    // > Edge
    int         *iso_edge_vtx  = NULL;
    PDM_g_num_t *iso_edge_gnum = NULL;
    int          iso_edge_n_group    = 0;
    int         *iso_edge_group_idx  = NULL;
    int         *iso_edge_group_lnum = NULL;
    PDM_g_num_t *iso_edge_group_gnum = NULL;
    int iso_n_edge = PDM_isosurface_connectivity_get(isos, id_iso, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                     NULL, &iso_edge_vtx,
                                                     PDM_OWNERSHIP_KEEP);
    PDM_isosurface_ln_to_gn_get(isos, id_iso, i_part, PDM_MESH_ENTITY_EDGE,
                               &iso_edge_gnum,
                                PDM_OWNERSHIP_KEEP);
    PDM_isosurface_group_get(isos, id_iso, i_part, PDM_MESH_ENTITY_EDGE,
                            &iso_edge_n_group, &iso_edge_group_idx, &iso_edge_group_lnum, &iso_edge_group_gnum,
                             PDM_OWNERSHIP_KEEP);
    // > Convert group into tag
    int         *iso_edge_tag1      = PDM_array_zeros_int (iso_n_edge);
    int         *iso_edge_tag2      = PDM_array_zeros_int (iso_n_edge);
    PDM_g_num_t *iso_edge_tag1_gnum = PDM_array_zeros_gnum(iso_n_edge);
    PDM_g_num_t *iso_edge_tag2_gnum = PDM_array_zeros_gnum(iso_n_edge);
    for (int i_group=0; i_group<iso_edge_n_group; ++i_group) {
      int i_beg_group = iso_edge_group_idx[i_group  ];
      int i_end_group = iso_edge_group_idx[i_group+1];
      for (int i_read=i_beg_group; i_read<i_end_group; ++i_read) {
        int         edge_lnum = iso_edge_group_lnum[i_read];
        PDM_g_num_t edge_gnum = iso_edge_group_gnum[i_read];
        if (iso_edge_tag1[edge_lnum-1]==0) {
          iso_edge_tag1     [edge_lnum-1] = i_group+1;
          iso_edge_tag1_gnum[edge_lnum-1] = edge_gnum;
        }
        iso_edge_tag2     [edge_lnum-1] = i_group+1;
        iso_edge_tag2_gnum[edge_lnum-1] = edge_gnum;
      }
    }


    // > Face
    int         *iso_face_vtx_idx  = NULL;
    int         *iso_face_vtx      = NULL;
    PDM_g_num_t *iso_face_gnum     = NULL;
    int iso_n_face = PDM_isosurface_connectivity_get(isos, id_iso, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX, &iso_face_vtx_idx, &iso_face_vtx, PDM_OWNERSHIP_KEEP);
    PDM_isosurface_ln_to_gn_get(isos, id_iso, i_part, PDM_MESH_ENTITY_FACE, &iso_face_gnum, PDM_OWNERSHIP_KEEP);


    /**
     * Write isosurface
     */
    char filename[999];

    const char *edge_field_name [] = {"tag1", "tag1_gnum", "tag2", "tag2_gnum"};
    const int  *edge_field_value[] = {iso_edge_tag1, iso_edge_tag1_gnum, iso_edge_tag2, iso_edge_tag2_gnum};
    sprintf(filename, "iso_edge_id_%i_part_%i_rank_%i.vtk", id_iso, i_part, i_rank);
    PDM_vtk_write_std_elements(filename,
                               iso_n_vtx,
                               iso_vtx_coord,
                               iso_vtx_gnum,
                               PDM_MESH_NODAL_BAR2,
                               iso_n_edge,
                               iso_edge_vtx,
                               iso_edge_gnum,
                               4,
                               edge_field_name,
                               edge_field_value);

    sprintf(filename, "iso_face_id_%i_part_%i_rank_%i.vtk", id_iso, i_part, i_rank);
    double *_iso_vtx_field = NULL;
    if (iso_vtx_field[id_iso]!=NULL) {
      _iso_vtx_field = iso_vtx_field[id_iso][i_part];
    }
    PDM_vtk_write_polydata_field(filename,
                                 iso_n_vtx,
                                 iso_vtx_coord,
                                 iso_vtx_gnum,
                                 iso_n_face,
                                 iso_face_vtx_idx,
                                 iso_face_vtx,
                                 iso_face_gnum,
                                 NULL,
                                 NULL,
                                 "field",
                                 _iso_vtx_field);

    free(iso_edge_tag1);
    free(iso_edge_tag2);
    free(iso_edge_tag1_gnum);
    free(iso_edge_tag2_gnum);
  }

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
     "  -nz      <level>  Number of vertices on the cube side (z direction).\n\n"
     "  -l       <level>  Cube length.\n\n"
     "  -n_part  <level>  Number of partitions (if partitioned entry).\n\n"
     "  -is_dist          Is entry distributed ou partitioned.\n\n"
     "  -local            Activate isosurface redistribution.\n\n"
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
  PDM_g_num_t          n_vtx_x    = 10;
  PDM_g_num_t          n_vtx_y    = 10;
  PDM_g_num_t          n_vtx_z    = 10;
  double               length     = 3.;
  int                  n_part     = 1;
  int                  order      = 1;
  PDM_Mesh_nodal_elt_t elt_type   = PDM_MESH_NODAL_TETRA4;
  int                  dist_entry = 0;
  int                  local      = 0;
  int                  visu       = 0;
  
  _read_args(argc,
             argv,
             &n_vtx_x,
             &n_vtx_y,
             &n_vtx_z,
             &length,
             &n_part,
             &dist_entry,
             &local,
             &visu);


  /*
   *  Generate mesh
   */
  PDM_dcube_nodal_t     *dcube_nodal = NULL;
  PDM_dmesh_nodal_t     *dmn         = NULL;
  PDM_part_mesh_nodal_t *pmn         = NULL;

  double  *iso_dfield = NULL;
  double **iso_field  = NULL;
  double **itp_field  = NULL;
  
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
                                      -0.5*length,
                                      -0.5*length,
                                      -0.5*length,
                                      length,
                                      length,
                                      length,
                                      n_vtx_x,
                                      n_vtx_y,
                                      n_vtx_z,
                                      n_part,
                                      PDM_SPLIT_DUAL_WITH_PARMETIS); // TODO: Allow various partitioning ?
    
    // > Fields initialisation
    iso_field = malloc(sizeof(double *) * n_part);
    itp_field = malloc(sizeof(double *) * n_part);

    for (int i_part=0; i_part<n_part; ++i_part) {
      int     n_vtx     = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);
      iso_field[i_part] = malloc(sizeof(double) * n_vtx);
      itp_field[i_part] = malloc(sizeof(double) * n_vtx);
      
      for (int i_vtx=0; i_vtx<n_vtx; ++i_vtx) {
        double c_4d[4] = {-0.08, 0.0, -0.8, -0.03};
        for (int i = 0; i < 4; i++) { c_4d[i] += 0.2; }
        iso_field[i_part][i_vtx] = _eval_julia4d(vtx_coord[3*i_vtx  ],
                                                 vtx_coord[3*i_vtx+1],
                                                 vtx_coord[3*i_vtx+2],
                                                 c_4d,
                                                 0);
        double c_2d[4] = {0.285, 0.01, 0., 0.};
        itp_field[i_part][i_vtx] = _eval_julia4d(vtx_coord[3*i_vtx  ],
                                                 vtx_coord[3*i_vtx+1],
                                                 vtx_coord[3*i_vtx+2],
                                                 c_2d,
                                                 1);
      }
    }
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

  // > User field isosurface

  double field_isovalues[1] = {0.1};
  int iso2 = PDM_isosurface_add(isos, 
                                PDM_ISO_SURFACE_KIND_FIELD,
                                1,
                                field_isovalues);
  if (dist_entry==1) {
    PDM_isosurface_dfield_set(isos,
                              iso2,
                              iso_dfield);
  }
  else if (dist_entry==0) {
    for (int i_part=0; i_part<n_part; ++i_part) {
      PDM_isosurface_field_set(isos,
                               iso2,
                               i_part,
                               iso_field[i_part]);
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_t_isosurface_2d_nodal dist_entry must be 0 or 1 (here set to %d)\n", dist_entry);
  }

  /*
   *  Compute isosurface
   */
  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_reset(isos, iso1);
  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_compute(isos, iso2);
  int n_iso = iso2+1;


  /*
   *  Interpolate field
   */
  double ***iso_itp_field = malloc(n_iso * sizeof(double **));
  for (int i_iso=0; i_iso<n_iso; ++i_iso) {
    iso_itp_field[i_iso] = NULL;
  } 
  if (dist_entry==1) {
    PDM_error(__FILE__, __LINE__, 0, "Dist entry not implemented\n");
  }
  else if (dist_entry==0) {
    if (local==1) {
      iso_itp_field[iso2] = malloc(n_part * sizeof(double **));
      for (int i_part=0; i_part<n_part; ++i_part) {
        
        int    *vtx_parent_idx  = NULL;
        int    *vtx_parent_lnum = NULL;
        double *vtx_parent_wght = NULL;
        int iso_n_vtx = PDM_isosurface_local_parent_get(isos, iso2, i_part, PDM_MESH_ENTITY_VTX, 
                                                       &vtx_parent_idx, &vtx_parent_lnum,
                                                        PDM_OWNERSHIP_KEEP);
        PDM_isosurface_vtx_parent_weight_get(isos, iso2, i_part, &vtx_parent_wght, PDM_OWNERSHIP_KEEP);
        
        iso_itp_field[iso2][i_part] = PDM_array_zeros_double(iso_n_vtx);
        for (int i_iso_vtx=0; i_iso_vtx<iso_n_vtx; ++i_iso_vtx) {
          int i_beg_parent = vtx_parent_idx[i_iso_vtx  ];
          int i_end_parent = vtx_parent_idx[i_iso_vtx+1];
          for (int i_parent=i_beg_parent; i_parent<i_end_parent; ++i_parent) {
            int    parent_lnum = vtx_parent_lnum[i_parent];
            double parent_wght = vtx_parent_wght[i_parent];
            iso_itp_field[iso2][i_part][i_iso_vtx] += itp_field[i_part][parent_lnum-1]*parent_wght;
          }
        }
      }
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Part entry with local=0 not implemented\n");
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_t_isosurface_2d_nodal dist_entry must be 0 or 1 (here set to %d)\n", dist_entry);
  }


  /*
   *  Visu isosurfaces
   */
  if (visu==1) {
    if (dist_entry==1) {
      PDM_error(__FILE__, __LINE__, 0, "PDM_t_isosurface_3d_nodal:: Not implemented\n");
    }
    else if (dist_entry==0) {
      // > iso line output
      PDM_part_mesh_nodal_t *iso_pmn1 = NULL;
      _build_pmn_from_iso_result(isos, iso1, n_part, &iso_pmn1, comm);

      PDM_part_mesh_nodal_dump_vtk(iso_pmn1,
                                   PDM_GEOMETRY_KIND_RIDGE,
                                   "pmn_iso_line_iso_mesh");
      PDM_part_mesh_nodal_dump_vtk(iso_pmn1,
                                   PDM_GEOMETRY_KIND_SURFACIC,
                                   "pmn_iso_face_iso_mesh");

      PDM_part_mesh_nodal_free(iso_pmn1);
      // _output_iso_result(isos, iso1, n_part, comm);
      _output_iso_result(isos, iso2, n_part, iso_itp_field, comm);
    }
  }


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

  if (iso_field!=NULL) {
    for (int i_part=0; i_part<n_part; ++i_part) {
      if (iso_field[i_part]!=NULL) {
        free(iso_field[i_part]);
      }
    }
    free(iso_field);
  }
  if (itp_field!=NULL) {
    for (int i_part=0; i_part<n_part; ++i_part) {
      if (itp_field[i_part]!=NULL) {
        free(itp_field[i_part]);
      }
    }
    free(itp_field);
  }
  for (int id_iso=0; id_iso<n_iso; ++id_iso) {
    if (iso_itp_field[id_iso]!=NULL) {
      for (int i_part=0; i_part<n_part; ++i_part) {
        if (iso_itp_field[id_iso][i_part]!=NULL) {
          free(iso_itp_field[id_iso][i_part]);
        }
      }
      free(iso_itp_field[id_iso]);
    }
  }
  free(iso_itp_field);


  PDM_MPI_Finalize();

  return 0;
}
