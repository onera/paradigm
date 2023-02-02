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
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_mesh_intersection.h"
#include "pdm_multipart.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_sphere_surf_gen.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

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
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -nA     <level>  Number vtx in side of mesh A (default : 10).\n\n"
     "  -nA     <level>  Number vtx in side of mesh B (default : 10).\n\n"
     "  -n_part <level>  Number vtx in side of mesh B (default : 10).\n\n"
     "  -t               Element kind .\n\n"
     "  -h               This message.\n\n");
  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int                    argc,
 char                 **argv,
 PDM_g_num_t           *n_vtx_a,
 PDM_g_num_t           *n_vtx_b,
 int                   *n_part,
 PDM_Mesh_nodal_elt_t  *elt_type
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-nA") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_vtx_a = atol(argv[i]);
        *n_vtx_a = (PDM_g_num_t) _n_vtx_a;
      }
    }
    else if (strcmp(argv[i], "-nB") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_vtx_b = atol(argv[i]);
        *n_vtx_b = (PDM_g_num_t) _n_vtx_b;
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
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}


static void _rotate_coord(const double angle,
                                double *coord)
{
  double Rz[3][3] = {{cos(angle), -sin(angle), 0},
                     {sin(angle),  cos(angle), 0},
                     {0            ,  0            , 1}};

  double x = coord[0];
  double y = coord[1];
  double z = coord[2];

  for (int j = 0; j < 3; j++) {
    coord[j] = Rz[j][0]*x + Rz[j][1]*y + Rz[j][2]*z;
  }
}



static void
_add_depth
(
 double *coord
 )
{
  double x = coord[0];
  double y = coord[1];
  double z = coord[2];

  // coord[1] = 0.8*y - 0.6*z;
  // coord[2] = 0.6*y + 0.8*z;

  // angular sector
  double t = PDM_PI * (0.5 + (x - 0.5 + 0.3*cos(3*y)) / 6.);
  double r = 0.3 + 0.6 * y;
  double scale = 0.2;
  // x = r * cos(t);
  // y = r * sin(t);
  coord[0] = 0.07*y*(1-y)*cos(2*PDM_PI*x) + 0.05*sin(5*y);//;0.5 * scale * (cos(3*(x+y) + .2) + sin(5*y + .1));
  coord[1] = r * cos(t);
  coord[2] = r * sin(t);

}


static
void
_generate_surface_mesh
(
 const PDM_MPI_Comm           comm,
 const PDM_g_num_t            n_vtx_seg,
 const PDM_Mesh_nodal_elt_t   elt_type,
 const int                    rotate,
 const double                 xmin,
 const double                 ymin,
 const double                 zmin,
 const double                 length,
 const PDM_split_dual_t       part_method,
 const int                    n_part,
       PDM_dmesh_nodal_t    **_dmn,
       PDM_multipart_t      **_mpart
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_dmesh_nodal_t *dmn = NULL;
  if (0) {
    PDM_sphere_surf_icosphere_gen_nodal(comm,
                                        n_vtx_seg,
                                        0, 0, 0,
                                        1,
                                        &dmn);
  }
  else {
    PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create (comm,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           length,
                                                           xmin,
                                                           ymin,
                                                           zmin,
                                                           elt_type,
                                                           1,
                                                           PDM_OWNERSHIP_USER);
    PDM_dcube_nodal_gen_build (dcube);
    dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
    PDM_dmesh_nodal_generate_distribution(dmn);
    PDM_dcube_nodal_gen_free(dcube);
  }

  PDM_g_num_t* distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];
  double* vtx_coord = PDM_DMesh_nodal_vtx_get(dmn);

  // randomize
  if (1) {
    double noise = 0.2*length/(double) (n_vtx_seg - 1);
    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      if (PDM_ABS(vtx_coord[3*i_vtx  ] - xmin         ) > 1.e-9 &&
          PDM_ABS(vtx_coord[3*i_vtx  ] - xmin - length) > 1.e-9 &&
          PDM_ABS(vtx_coord[3*i_vtx+1] - ymin         ) > 1.e-9 &&
          PDM_ABS(vtx_coord[3*i_vtx+1] - ymin - length) > 1.e-9) {
        srand(distrib_vtx[i_rank] + i_vtx);
        for (int i = 0; i < 2; i++) {
          vtx_coord[3*i_vtx+i] += noise*0.5*(2*rand()/(double) RAND_MAX - 1);
        }
      }
    }
  }

  if(rotate) {
    // Do something
    double pi = 4 * atan(1.);
    double angle = pi/5.;
    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      _rotate_coord(angle, &vtx_coord[3*i_vtx]);
    }
  }

  if (1) {
    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      _add_depth(&vtx_coord[3*i_vtx]);
    }
  }

  if(0 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "sphere_surf_");
  }

  int n_zone = 1;
  // int n_part_zones = {n_part};
  int *n_part_zones = (int *) malloc(sizeof(int) * n_zone);
  n_part_zones[0] = n_part;

  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                n_part_zones,
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

  PDM_multipart_register_dmesh_nodal(mpart, 0, dmn);
  PDM_multipart_run_ppart(mpart);

  free(n_part_zones);


  *_mpart = mpart;
  *_dmn   = dmn;

}

static
void
_set_mesh
(
 PDM_mesh_intersection_t *mi,
 PDM_ol_mesh_t            i_mesh,
 PDM_multipart_t         *mpart,
 int                      n_part
)
{
  for (int i_part = 0; i_part < n_part; i_part++) {

    int *face_edge_idx;
    int *face_edge;
    int *edge_vtx_idx;
    int *edge_vtx;

    int n_proc, tn_part;
    int _n_vtx, n_bounds, n_joins, n_part_joins;
    int sface_edge, sedge_vtx, sedge_bound, sedge_join;
    int  n_section;
    int* n_elt;

    int n_face, n_edge;
    PDM_multipart_part_dim_get(mpart, 0, i_part, &n_section, &n_elt,
                               &n_face, &n_edge, &n_part_joins, &_n_vtx, &n_proc, &tn_part,
                               &sface_edge, &sedge_vtx, &sedge_bound, &n_bounds, &sedge_join, &n_joins);

    double       *_vtx;
    int          *_edge_face;
    int          *edge_bound_idx, *edge_bound, *edge_join_idx, *edge_join;
    int          *edge_part_bound_proc_idx, *edge_part_bound_part_idx, *edge_part_bound;
    PDM_g_num_t  *_face_ln_to_gn, *edge_ln_to_gn, *_vtx_ln_to_gn, *edge_bound_ln_to_gn, *edge_join_ln_to_gn;
    int          *face_tag, *edge_tag, *vtx_tag;
    int         **elt_vtx_idx;
    int         **elt_vtx;
    PDM_g_num_t **elt_section_ln_to_gn;

    PDM_multipart_part_val_get(mpart, 0, i_part, &elt_vtx_idx, &elt_vtx, &elt_section_ln_to_gn,
                               &face_tag, &face_edge_idx, &face_edge, &_face_ln_to_gn,
                               &edge_tag, &_edge_face, &edge_vtx_idx, &edge_vtx, &edge_ln_to_gn,
                               &edge_part_bound_proc_idx, &edge_part_bound_part_idx, &edge_part_bound,
                               &vtx_tag, &_vtx, &_vtx_ln_to_gn, &edge_bound_idx, &edge_bound,
                               &edge_bound_ln_to_gn, &edge_join_idx, &edge_join, &edge_join_ln_to_gn);

    double *vtx_coord;
    int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                 0,
                                                 i_part,
                                                 &vtx_coord,
                                                 PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *face_ln_to_gn;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    &face_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_EDGE,
                                    &edge_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);
    PDM_g_num_t *vtx_ln_to_gn;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_VERTEX,
                                    &vtx_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    n_face = PDM_multipart_part_connectivity_get(mpart,
                                                 0,
                                                 i_part,
                                                 PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                 &face_edge,
                                                 &face_edge_idx,
                                                 PDM_OWNERSHIP_KEEP);
    n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                 0,
                                                 i_part,
                                                 PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                 &edge_vtx,
                                                 &edge_vtx_idx,
                                                 PDM_OWNERSHIP_KEEP);

    PDM_mesh_intersection_part_set(mi,
                                   i_mesh,
                                   i_part,
                                   0, //n_cell
                                   n_face,
                                   n_edge,
                                   n_vtx,
                                   NULL,//cell_face_idx,
                                   NULL,//cell_face,
                                   face_edge_idx,
                                   face_edge,
                                   edge_vtx,
                                   NULL,//face_vtx_idx,
                                   NULL,//face_vtx,
                                   NULL,//cell_ln_to_gn,
                                   face_ln_to_gn,
                                   edge_ln_to_gn,
                                   vtx_ln_to_gn,
                                   vtx_coord);
  }

}


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  PDM_g_num_t n_vtx_a   = 10;
  PDM_g_num_t n_vtx_b   = 10;
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_TRIA3;

#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  int n_part = 1;

  _read_args(argc,
             argv,
             &n_vtx_a,
             &n_vtx_b,
             &n_part,
             &elt_type);

  /*
   * Generate meshA
   */
  double length_a = 1.;
  int rotate_a = 0;
  PDM_dmesh_nodal_t     *dmn_surf_a   = NULL;
  PDM_multipart_t       *mpart_surf_a = NULL;
  _generate_surface_mesh (comm,
                          n_vtx_a,
                          elt_type,
                          rotate_a,
                          0.,
                          0.,
                          0.,
                          length_a,
                          part_method,
                          n_part,
                          &dmn_surf_a,
                          &mpart_surf_a);


  double length_b = 1.;
  int rotate_b = 1;
  PDM_dmesh_nodal_t     *dmn_surf_b   = NULL;
  PDM_multipart_t       *mpart_surf_b = NULL;
  _generate_surface_mesh (comm,
                          n_vtx_b,
                          elt_type,
                          rotate_b,
                          0.5,
                          0.5,
                          0.,
                          length_b,
                          part_method,
                          n_part,
                          &dmn_surf_b,
                          &mpart_surf_b);

  if(1 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn_surf_a,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "dmn_surf_a_");
    PDM_dmesh_nodal_dump_vtk(dmn_surf_b,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "dmn_surf_b_");
  }

  /*
   * Mesh_intersection
   */
  int dim_mesh_a = 2;
  int dim_mesh_b = 2;
  PDM_mesh_intersection_t* mi = PDM_mesh_intersection_create(PDM_MESH_INTERSECTION_KIND_SOFT,
                                                             dim_mesh_a,
                                                             dim_mesh_b,
                                                             n_part,
                                                             n_part,
                                                             1e-6,
                                                             comm);

  /*
   * Set mesh_a and mesh_b
   */
  _set_mesh(mi, PDM_OL_MESH_A, mpart_surf_a, n_part);
  _set_mesh(mi, PDM_OL_MESH_B, mpart_surf_b, n_part);

  PDM_mesh_intersection_compute(mi);

  PDM_mesh_intersection_free(mi);


  PDM_DMesh_nodal_free(dmn_surf_b);
  PDM_multipart_free(mpart_surf_b);

  PDM_DMesh_nodal_free(dmn_surf_a);
  PDM_multipart_free(mpart_surf_a);

  PDM_MPI_Barrier(comm);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize ();

  return 0;

}
