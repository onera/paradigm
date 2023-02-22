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
#include "pdm_gnum.h"
#include "pdm_part_connectivity_transform.h"

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
     "  -n      <level>  Number vtx in side of mesh A (default : 10).\n\n"
     "  -n_part <level>  Number of partition                        .\n\n"
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


static
void
_generate_volume_mesh
(
 const PDM_MPI_Comm           comm,
 const PDM_g_num_t            n_vtx_seg,
 const PDM_Mesh_nodal_elt_t   elt_type,
 const int                    rotate,
 const double                 xmin,
 const double                 ymin,
 const double                 zmin,
 const double                 lenght,
 const PDM_split_dual_t       part_method,
 const int                    n_part,
       PDM_dmesh_nodal_t    **_dmn,
       PDM_multipart_t      **_mpart
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         lenght,
                                                         xmin,
                                                         ymin,
                                                         zmin,
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_USER);
  PDM_dcube_nodal_gen_build (dcube);
  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);
  PDM_dcube_nodal_gen_free(dcube);

  PDM_UNUSED(rotate);
  // if(rotate) {
  //   // Do something
  //   double pi = 4 * atan(1.);
  //   double angle = pi/5.;
  //   PDM_g_num_t* distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  //   int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];
  //   double* vtx_coord = PDM_DMesh_nodal_vtx_get(dmn);
  //   for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
  //     _rotate_coord(angle, &vtx_coord[3*i_vtx]);
  //   }
  // }

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
_create_wall_surf
(
 const PDM_MPI_Comm           comm,
 const int                    n_part,
       PDM_multipart_t       *mpart,
       double              ***psurf_vtx_coord_out,
       int                 ***psurf_face_vtx_idx_out,
       int                 ***psurf_face_vtx_out,
       PDM_g_num_t         ***psurf_face_ln_to_gn_out,
       PDM_g_num_t         ***psurf_vtx_ln_to_gn_out
)
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  int          *n_surf_vtx           = malloc(n_part * sizeof(int          ));
  int          *n_surf_face          = malloc(n_part * sizeof(int          ));

  double      **psurf_vtx_coord      = malloc(n_part * sizeof(double      *));
  int         **psurf_face_vtx_idx   = malloc(n_part * sizeof(int         *));
  int         **psurf_face_vtx       = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **psurf_face_ln_to_gn  = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **psurf_vtx_ln_to_gn   = malloc(n_part * sizeof(PDM_g_num_t *));

  /* Compute gnum for vtx and faces */
  PDM_gen_gnum_t* gnum_face = PDM_gnum_create(3,
                                              n_part,
                                              PDM_FALSE,
                                              1.e-6,
                                              comm,
                                              PDM_OWNERSHIP_USER);

  PDM_gen_gnum_t* gnum_vtx = PDM_gnum_create(3,
                                             n_part,
                                             PDM_FALSE,
                                             1.e-6,
                                             comm,
                                             PDM_OWNERSHIP_USER);

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

    int  n_bound = 0;
    int* group_face_idx      = NULL;
    int* group_face          = NULL;
    PDM_g_num_t* face_group_ln_to_gn = NULL;

    PDM_multipart_bound_get(mpart,
                            0,
                            i_part,
                            PDM_BOUND_TYPE_FACE,
                            &n_bound,
                            &group_face_idx,
                            &group_face,
                            &face_group_ln_to_gn);

    int *face_vtx = NULL;
    int *face_vtx_idx = face_edge_idx;
    PDM_compute_face_vtx_from_face_and_edge(n_face, face_edge_idx, face_edge, edge_vtx, &face_vtx);

    /*
     * Nez de la plaque en x = 0
     */
    int i_group = 0;
    n_surf_face[i_part] = 0;
    int n_surf_face_vtx = 0;
    int n_face_in_group = group_face_idx[i_group+1] - group_face_idx[i_group];
    int* face_bnd = malloc(n_face_in_group * sizeof(int));
    for(int idx_face = group_face_idx[i_group]; idx_face < group_face_idx[i_group+1]; ++idx_face) {
      int i_face = group_face[idx_face]-1;

      int is_in_plate = 1;
      for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = face_vtx[idx_vtx]-1;
        if(vtx_coord[3*i_vtx] < 0.) {
          is_in_plate = 0;
        }
      }

      if(is_in_plate == 0) {
        continue;
      }

      face_bnd[n_surf_face[i_part]++] = i_face;

      n_surf_face_vtx += face_vtx_idx[i_face+1] - face_vtx_idx[i_face];

    }


    psurf_vtx_coord    [i_part] = malloc(3 * n_surf_face_vtx          * sizeof(double     ));
    psurf_face_vtx_idx [i_part] = malloc(   ( n_surf_face[i_part] +1) * sizeof(int        ));
    psurf_face_vtx     [i_part] = malloc(    n_surf_face_vtx          * sizeof(int        ));
    psurf_face_ln_to_gn[i_part] = malloc(    n_surf_face[i_part]      * sizeof(PDM_g_num_t));
    psurf_vtx_ln_to_gn [i_part] = malloc(    n_surf_face_vtx          * sizeof(PDM_g_num_t));

    double      *_psurf_vtx_coord     = psurf_vtx_coord    [i_part];
    int         *_psurf_face_vtx_idx  = psurf_face_vtx_idx [i_part];
    int         *_psurf_face_vtx      = psurf_face_vtx     [i_part];
    PDM_g_num_t *_psurf_face_ln_to_gn = psurf_face_ln_to_gn[i_part];
    PDM_g_num_t *_psurf_vtx_ln_to_gn  = psurf_vtx_ln_to_gn [i_part];

    int *vtx_flags = malloc(n_vtx * sizeof(int));
    for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx) {
      vtx_flags[i_vtx] = -100;
    }
    n_surf_vtx[i_part] = 0;

    _psurf_face_vtx_idx[0] = 0;
    for(int idx_face = 0; idx_face < n_surf_face[i_part]; ++idx_face) {
      int i_face = face_bnd[idx_face];

      _psurf_face_vtx_idx[idx_face+1] = _psurf_face_vtx_idx[idx_face];
      _psurf_face_ln_to_gn[idx_face] = face_ln_to_gn[i_face];

      for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = face_vtx[idx_vtx]-1;

        if(vtx_flags[i_vtx] == -100) {
          vtx_flags[i_vtx] = n_surf_vtx[i_part];
          _psurf_vtx_coord[3*n_surf_vtx[i_part]  ] = vtx_coord[3*i_vtx  ];
          _psurf_vtx_coord[3*n_surf_vtx[i_part]+1] = vtx_coord[3*i_vtx+1];
          _psurf_vtx_coord[3*n_surf_vtx[i_part]+2] = vtx_coord[3*i_vtx+2];

          _psurf_vtx_ln_to_gn[n_surf_vtx[i_part]] = vtx_ln_to_gn[i_vtx];

          n_surf_vtx[i_part]++;
        }

        _psurf_face_vtx[_psurf_face_vtx_idx[idx_face+1]++] = vtx_flags[i_vtx]+1;
      }

    }

    printf("n_surf_face = %i \n", n_surf_face[i_part]);
    printf("n_surf_vtx  = %i \n", n_surf_vtx [i_part]);

    psurf_vtx_coord    [i_part] = realloc(psurf_vtx_coord    [i_part], 3 * n_surf_vtx[i_part] * sizeof(double     ));
    psurf_vtx_ln_to_gn [i_part] = realloc(psurf_vtx_ln_to_gn [i_part],     n_surf_vtx[i_part] * sizeof(PDM_g_num_t));

    PDM_gnum_set_from_parents(gnum_face, i_part, n_surf_face[i_part], _psurf_face_ln_to_gn);
    PDM_gnum_set_from_parents(gnum_vtx , i_part, n_surf_vtx [i_part], _psurf_vtx_ln_to_gn );

    free(vtx_flags);
    free(face_bnd);
    free(face_vtx);

  }

  PDM_gnum_compute(gnum_face);
  PDM_gnum_compute(gnum_vtx );

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(psurf_face_ln_to_gn[i_part]);
    free(psurf_vtx_ln_to_gn [i_part]);
    psurf_face_ln_to_gn[i_part] = PDM_gnum_get(gnum_face, 0);
    psurf_vtx_ln_to_gn [i_part] = PDM_gnum_get(gnum_vtx , 0);
  }

  PDM_gnum_free(gnum_face);
  PDM_gnum_free(gnum_vtx);


  /* Vtk en légende */
  if(1 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      char filename[999];
      sprintf(filename, "face_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_polydata(filename,
                             n_surf_vtx [i_part],
                             psurf_vtx_coord    [i_part],
                             psurf_vtx_ln_to_gn [i_part],
                             n_surf_face        [i_part],
                             psurf_face_vtx_idx [i_part],
                             psurf_face_vtx     [i_part],
                             psurf_face_ln_to_gn[i_part],
                             NULL);
    }

  }




  *psurf_vtx_coord_out      = psurf_vtx_coord;
  *psurf_face_vtx_idx_out   = psurf_face_vtx_idx;
  *psurf_face_vtx_out       = psurf_face_vtx;
  *psurf_face_ln_to_gn_out  = psurf_face_ln_to_gn;
  *psurf_vtx_ln_to_gn_out   = psurf_vtx_ln_to_gn;
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
  PDM_MPI_Comm_rank (comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_g_num_t n_vtx_a   = 10;
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_HEXA8;

  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_HILBERT;

  int n_part = 1;

  _read_args(argc,
             argv,
             &n_vtx_a,
             &n_part,
             &elt_type);


  /*
   * Generate meshA
   */
  double lenght_a = 1.;
  int rotate_a = 0;
  PDM_dmesh_nodal_t     *dmn_vol_a   = NULL;
  PDM_multipart_t       *mpart_vol_a = NULL;
  _generate_volume_mesh (comm,
                         n_vtx_a,
                         elt_type,
                         rotate_a,
                         0.,
                         0.,
                         0.,
                         lenght_a,
                         part_method,
                         n_part,
                         &dmn_vol_a,
                         &mpart_vol_a);

  if(1 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn_vol_a,
                             PDM_GEOMETRY_KIND_VOLUMIC,
                             "dmn_vol_a_");
    PDM_dmesh_nodal_dump_vtk(dmn_vol_a,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "dmn_surf_a_");

  }

  /*
   * Extract boundary wall
   */
  double      **psurf_vtx_coord      = NULL;
  int         **psurf_face_vtx_idx   = NULL;
  int         **psurf_face_vtx       = NULL;
  PDM_g_num_t **psurf_face_ln_to_gn  = NULL;
  PDM_g_num_t **psurf_vtx_ln_to_gn   = NULL;
  _create_wall_surf(comm,
                    n_part,
                    mpart_vol_a,
                    &psurf_vtx_coord,
                    &psurf_face_vtx_idx,
                    &psurf_face_vtx,
                    &psurf_face_ln_to_gn,
                    &psurf_vtx_ln_to_gn );



  PDM_DMesh_nodal_free(dmn_vol_a);
  PDM_multipart_free(mpart_vol_a);
  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize ();

  return 0;

}
