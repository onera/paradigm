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
#include "pdm_dist_cloud_surf.h"
#include "pdm_gnum.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dbbtree.h"
#include "pdm_geom_elem.h"

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
                                                         6,
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


  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

  double pi = 4 * atan(1.);
  for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {

    double x = dvtx_coord[3*i_vtx  ];
    double y = dvtx_coord[3*i_vtx+1];
    double z = dvtx_coord[3*i_vtx+2];

    dvtx_coord[3*i_vtx+1] = tanh( y * y * y);
    y = dvtx_coord[3*i_vtx+1];

    double angle = -pi/4;
    double Rz[3][3] = {{cos(angle), -sin(angle), 0},
                       {sin(angle),  cos(angle), 0},
                       {0         ,  0         , 1}};

    // dvtx_coord[3*i_vtx  ] = x +
    // dvtx_coord[3*i_vtx+1] = cos(y) - sin(x);

    if( x > 0.) {
      dvtx_coord[3*i_vtx  ] *= 4;
    }
    x = dvtx_coord[3*i_vtx];

    for (int j = 0; j < 1; j++) {
      dvtx_coord[3*i_vtx+j] = Rz[j][0]*x + Rz[j][1]*y + Rz[j][2]*z;
    }


  }

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
_cell_center_3d
(
  int      pn_cell,
  int     *pcell_face_idx,
  int     *pcell_face,
  int     *pface_edge_idx,
  int     *pface_edge,
  int     *pface_vtx_idx,
  int     *pface_vtx,
  int     *pedge_vtx,
  double  *pvtx_coord,
  double **cell_center
)
{
  int from_edge = 0;
  int from_face = 0;
  if(pface_edge     != NULL) {
    from_edge = 1;
  }
  if(pface_vtx     != NULL) {
    from_face = 1;
  }
  assert(pvtx_coord     != NULL);

  double* entity_center = malloc(3 * pn_cell * sizeof(double ));

  if(from_face == 1) {
    for(int i_cell = 0; i_cell < pn_cell; ++i_cell) {

      entity_center[3*i_cell  ] = 0.;
      entity_center[3*i_cell+1] = 0.;
      entity_center[3*i_cell+2] = 0.;

      double inv = 1./((double) pcell_face_idx[i_cell+1] - pcell_face_idx[i_cell]);

      for(int idx_face = pcell_face_idx[i_cell]; idx_face < pcell_face_idx[i_cell+1]; ++idx_face) {
        int i_face = PDM_ABS(pcell_face[idx_face])-1;

        double inv2 = 1./((double)  pface_vtx_idx[i_face+1] - pface_vtx_idx[i_face]);

        double fcx = 0;
        double fcy = 0;
        double fcz = 0;
        for(int idx_vtx = pface_vtx_idx[i_face]; idx_vtx < pface_vtx_idx[i_face+1]; ++idx_vtx) {
          int i_vtx = pface_vtx[idx_vtx]-1;
          fcx += pvtx_coord[3*i_vtx  ];
          fcy += pvtx_coord[3*i_vtx+1];
          fcz += pvtx_coord[3*i_vtx+2];
        }
        fcx = fcx * inv2;
        fcy = fcy * inv2;
        fcz = fcz * inv2;

        entity_center[3*i_cell  ] += fcx;
        entity_center[3*i_cell+1] += fcy;
        entity_center[3*i_cell+2] += fcz;
      }

      entity_center[3*i_cell  ] = entity_center[3*i_cell  ] * inv;
      entity_center[3*i_cell+1] = entity_center[3*i_cell+1] * inv;
      entity_center[3*i_cell+2] = entity_center[3*i_cell+2] * inv;
    } /* End cell */
  } else if( from_edge == 1) {
    for(int i_cell = 0; i_cell < pn_cell; ++i_cell) {

      entity_center[3*i_cell  ] = 0.;
      entity_center[3*i_cell+1] = 0.;
      entity_center[3*i_cell+2] = 0.;

      double inv = 1./((double)  pcell_face_idx[i_cell+1] - pcell_face_idx[i_cell]);

      double fcx = 0;
      double fcy = 0;
      double fcz = 0;
      for(int idx_face = pcell_face_idx[i_cell]; idx_face < pcell_face_idx[i_cell+1]; ++idx_face) {
        int i_face = PDM_ABS(pcell_face[idx_face])-1;

        double inv2 = 1./((double)  pface_edge_idx[i_face+1] - pface_edge_idx[i_face]);

        for(int idx_edge = pface_edge_idx[i_face]; idx_edge < pface_edge_idx[i_face+1]; ++idx_edge) {
          int i_edge = PDM_ABS(pface_edge[idx_edge])-1;
          int i_vtx1 = pedge_vtx[2*i_edge  ] - 1;
          int i_vtx2 = pedge_vtx[2*i_edge+1] - 1;
          fcx += 0.5 * (pvtx_coord[3*i_vtx1  ] + pvtx_coord[3*i_vtx2  ]);
          fcy += 0.5 * (pvtx_coord[3*i_vtx1+1] + pvtx_coord[3*i_vtx2+1]);
          fcz += 0.5 * (pvtx_coord[3*i_vtx1+2] + pvtx_coord[3*i_vtx2+2]);
        }
        fcx = fcx * inv2;
        fcy = fcy * inv2;
        fcz = fcz * inv2;

        entity_center[3*i_cell  ] += fcx;
        entity_center[3*i_cell+1] += fcy;
        entity_center[3*i_cell+2] += fcz;
      }

      entity_center[3*i_cell  ] = entity_center[3*i_cell  ] * inv;
      entity_center[3*i_cell+1] = entity_center[3*i_cell+1] * inv;
      entity_center[3*i_cell+2] = entity_center[3*i_cell+2] * inv;
    } /* End cell */
  }

  *cell_center = entity_center;
}

static
void
_create_wall_surf
(
 const PDM_MPI_Comm           comm,
 const int                    n_part,
       PDM_multipart_t       *mpart,
       int                  **n_surf_vtx_out,
       int                  **n_surf_face_out,
       double              ***psurf_vtx_coord_out,
       int                 ***psurf_face_vtx_idx_out,
       int                 ***psurf_face_vtx_out,
       PDM_g_num_t         ***psurf_face_ln_to_gn_out,
       PDM_g_num_t         ***psurf_vtx_ln_to_gn_out,
       int                  **pn_cell_out,
       PDM_g_num_t         ***pcell_ln_to_gn_out,
       double              ***cell_center_out
)
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  int          *n_surf_vtx           = malloc(n_part * sizeof(int          ));
  int          *n_surf_face          = malloc(n_part * sizeof(int          ));

  double      **cell_center          = malloc(n_part * sizeof(double      *));
  double      **psurf_vtx_coord      = malloc(n_part * sizeof(double      *));
  int         **psurf_face_vtx_idx   = malloc(n_part * sizeof(int         *));
  int         **psurf_face_vtx       = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **psurf_face_ln_to_gn  = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **psurf_vtx_ln_to_gn   = malloc(n_part * sizeof(PDM_g_num_t *));

  int          *pn_cell        = malloc(n_part * sizeof(int          ));
  PDM_g_num_t **pcell_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));

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

    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    &pcell_ln_to_gn[i_part],
                                    PDM_OWNERSHIP_KEEP);

    int *pcell_face     = NULL;
    int *pcell_face_idx = NULL;
    int n_cell =   PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &pcell_face,
                                                       &pcell_face_idx,
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

    pn_cell[i_part] = n_cell;
    _cell_center_3d(n_cell,
                    pcell_face_idx,
                    pcell_face,
                    NULL,
                    NULL,
                    face_vtx_idx,
                    face_vtx,
                    NULL,
                    vtx_coord,
                    &cell_center[i_part]);
    /*
     * Nez de la plaque en x = 0
     */
    int i_group = 4;
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
    PDM_gnum_set_from_parents(gnum_vtx , i_part, n_surf_vtx [i_part], psurf_vtx_ln_to_gn [i_part] );

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


  *pn_cell_out        = pn_cell;
  *pcell_ln_to_gn_out = pcell_ln_to_gn;
  *cell_center_out    = cell_center;
  *n_surf_vtx_out     = n_surf_vtx;
  *n_surf_face_out    = n_surf_face;

  *psurf_vtx_coord_out      = psurf_vtx_coord;
  *psurf_face_vtx_idx_out   = psurf_face_vtx_idx;
  *psurf_face_vtx_out       = psurf_face_vtx;
  *psurf_face_ln_to_gn_out  = psurf_face_ln_to_gn;
  *psurf_vtx_ln_to_gn_out   = psurf_vtx_ln_to_gn;
}


static
void
_create_wall_ray
(
 const PDM_MPI_Comm           comm,
 const int                    n_part,
       int                  *n_surf_vtx,
       int                  *n_surf_face,
       double              **psurf_vtx_coord,
       int                 **psurf_face_vtx_idx,
       int                 **psurf_face_vtx,
       PDM_g_num_t         **psurf_face_ln_to_gn,
       PDM_g_num_t         **psurf_vtx_ln_to_gn,
       int                  *pn_ray_out,
       PDM_g_num_t         **pray_ln_to_gn_out,
       double              **pray_coord_out
)
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  int pn_ray = 0;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_ray += n_surf_face[i_part];
  }

  PDM_g_num_t *pray_ln_to_gn = malloc(pn_ray * sizeof(PDM_g_num_t));
  double      *pray_coord    = malloc(6 * pn_ray * sizeof(double     ));

  pn_ray = 0;
  for(int i_part = 0; i_part < n_part; ++i_part) {


    double      *face_normal    = malloc(3 * n_surf_face[i_part] * sizeof(double     ));
    double      *face_center    = malloc(3 * n_surf_face[i_part] * sizeof(double     ));
    PDM_geom_elem_polygon_properties(n_surf_face[i_part],
                                     psurf_face_vtx_idx[i_part],
                                     psurf_face_vtx    [i_part],
                                     psurf_vtx_coord   [i_part],
                                     face_normal,
                                     face_center,
                                     NULL,
                                     NULL);

    double dmax = 0.5;
    for(int i_face = 0; i_face < n_surf_face[i_part]; ++i_face) {

      pray_coord[6*pn_ray  ] = face_center[3*i_face  ];
      pray_coord[6*pn_ray+1] = face_center[3*i_face+1];
      pray_coord[6*pn_ray+2] = face_center[3*i_face+2];

      double nx = face_normal[3*i_face  ];
      double ny = face_normal[3*i_face+1];
      double nz = face_normal[3*i_face+2];

      double sn = sqrt(nx * nx + ny * ny + nz * nz);
      nx = - nx / sn;
      ny = - ny / sn;
      nz = - nz / sn;

      double dnx = nx * dmax;
      double dny = ny * dmax;
      double dnz = nz * dmax;

      double xb = pray_coord[6*pn_ray  ] + dnx;
      double yb = pray_coord[6*pn_ray+1] + dny;
      double zb = pray_coord[6*pn_ray+2] + dnz;

      pray_coord[6*pn_ray+3] = xb;
      pray_coord[6*pn_ray+4] = yb;
      pray_coord[6*pn_ray+5] = zb;

      pray_ln_to_gn[pn_ray++] = psurf_face_ln_to_gn[i_face];

    }

    free(face_normal);
    free(face_center);
  }


  char filename[999];
  sprintf(filename, "ray_%i.vtk", i_rank);
  PDM_vtk_write_lines(filename,
                      pn_ray,
                      pray_coord,
                      pray_ln_to_gn,
                      NULL);


  *pn_ray_out        = pn_ray;
  *pray_ln_to_gn_out = pray_ln_to_gn;
  *pray_coord_out    = pray_coord;

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
                         -0.2,
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
  int          *pn_cell              = NULL;
  PDM_g_num_t **pcell_ln_to_gn       = NULL;
  int          *psurf_vtx            = NULL;
  int          *psurf_face           = NULL;
  double      **psurf_vtx_coord      = NULL;
  int         **psurf_face_vtx_idx   = NULL;
  int         **psurf_face_vtx       = NULL;
  PDM_g_num_t **psurf_face_ln_to_gn  = NULL;
  PDM_g_num_t **psurf_vtx_ln_to_gn   = NULL;
  double      **cell_center          = NULL;
  _create_wall_surf(comm,
                    n_part,
                    mpart_vol_a,
                    &psurf_vtx,
                    &psurf_face,
                    &psurf_vtx_coord,
                    &psurf_face_vtx_idx,
                    &psurf_face_vtx,
                    &psurf_face_ln_to_gn,
                    &psurf_vtx_ln_to_gn,
                    &pn_cell,
                    &pcell_ln_to_gn,
                    &cell_center);

  /* Wall distance */
  int n_point_cloud = 1;
  PDM_dist_cloud_surf_t* dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_MESH_SETTED,
                                                            n_point_cloud,
                                                            comm,
                                                            PDM_OWNERSHIP_KEEP);

  PDM_dist_cloud_surf_surf_mesh_global_data_set (dist,
                                                 n_part);

  PDM_dist_cloud_surf_n_part_cloud_set (dist, 0, n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    PDM_dist_cloud_surf_surf_mesh_part_set (dist,
                                            i_part,
                                            psurf_face         [i_part],
                                            psurf_face_vtx_idx [i_part],
                                            psurf_face_vtx     [i_part],
                                            psurf_face_ln_to_gn[i_part],
                                            psurf_vtx          [i_part],
                                            psurf_vtx_coord    [i_part],
                                            psurf_vtx_ln_to_gn [i_part]);

    PDM_dist_cloud_surf_cloud_set (dist,
                                   0,
                                   i_part,
                                   pn_cell       [i_part],
                                   cell_center   [i_part],
                                   pcell_ln_to_gn[i_part]);
  }

  PDM_dist_cloud_surf_compute(dist);

  for (int i_part = 0; i_part < n_part; i_part++) {

    double      *distance;
    double      *projected;
    PDM_g_num_t *closest_elt_gnum;

    PDM_dist_cloud_surf_get (dist,
                             0,
                             i_part,
                             &distance,
                             &projected,
                             &closest_elt_gnum);

    char filename[999];
    sprintf(filename, "distance_%3.3d_%3.3d.vtk", i_part, i_rank);

    const char   *vector_field_name[1] = {"distance"};
    const double *vector_field     [1] = {distance};
    PDM_vtk_write_point_cloud_with_field(filename,
                                         pn_cell       [i_part],
                                         cell_center   [i_part],
                                         pcell_ln_to_gn[i_part],
                                         NULL,
                                         1,
                                         vector_field_name,
                                         vector_field,
                                         0,
                                         NULL,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL );
  }

  /* Create field of speed */
  double **velocity = malloc(n_part * sizeof(double));
  for (int i_part = 0; i_part < n_part; i_part++) {

    double      *distance;
    double      *projected;
    PDM_g_num_t *closest_elt_gnum;

    PDM_dist_cloud_surf_get (dist,
                             0,
                             i_part,
                             &distance,
                             &projected,
                             &closest_elt_gnum);

    velocity[i_part] = malloc(pn_cell[i_part] * sizeof(double));

    for(int i_cell = 0; i_cell < pn_cell[i_part]; ++i_cell) {
      if(cell_center   [i_part][3*i_cell+1] < 0.1 * cell_center   [i_part][3*i_cell]) {
        velocity[i_part][i_cell] = tanh(6*distance[i_cell]);
      } else {
        velocity[i_part][i_cell] = 1.;
      }
    }

    char filename[999];
    sprintf(filename, "velocity_%3.3d_%3.3d.vtk", i_part, i_rank);

    const char   *vector_field_name[1] = {"velocity"};
    const double *vector_field     [1] = {velocity[i_part]};
    PDM_vtk_write_point_cloud_with_field(filename,
                                         pn_cell       [i_part],
                                         cell_center   [i_part],
                                         pcell_ln_to_gn[i_part],
                                         NULL,
                                         1,
                                         vector_field_name,
                                         vector_field,
                                         0,
                                         NULL,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL );
  }

  /*
   * Compute extents of all cell
   */
  double **box_extents = malloc(n_part * sizeof(double *));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int *pcell_face     = NULL;
    int *pcell_face_idx = NULL;
    int *pface_edge      = NULL;
    int *pface_edge_idx  = NULL;
    int *pedge_vtx      = NULL;
    int *pedge_vtx_idx  = NULL;

    double *vtx_coord;
    int n_vtx = PDM_multipart_part_vtx_coord_get(mpart_vol_a,
                                                 0,
                                                 i_part,
                                                 &vtx_coord,
                                                 PDM_OWNERSHIP_KEEP);

    int n_edge = PDM_multipart_part_connectivity_get(mpart_vol_a,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                     &pedge_vtx,
                                                     &pedge_vtx_idx,
                                                     PDM_OWNERSHIP_KEEP);

    int n_face = PDM_multipart_part_connectivity_get(mpart_vol_a,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                     &pface_edge,
                                                     &pface_edge_idx,
                                                     PDM_OWNERSHIP_KEEP);

    int n_cell =   PDM_multipart_part_connectivity_get(mpart_vol_a,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &pcell_face,
                                                       &pcell_face_idx,
                                                       PDM_OWNERSHIP_KEEP);

    int *pface_vtx     = NULL;
    int *pface_vtx_idx = pface_edge_idx;
    PDM_compute_face_vtx_from_face_and_edge(n_face,
                                            pface_edge_idx,
                                            pface_edge,
                                            pedge_vtx,
                                            &pface_vtx);

    box_extents[i_part] = malloc(6 * n_cell * sizeof(double));

    for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
      double *_extents = box_extents[i_part] + 6 * i_cell;
      for(int k = 0; k < 3; ++k) {
        _extents[k  ] =  HUGE_VAL;
        _extents[k+3] = -HUGE_VAL;
      }
      for(int idx_face = pcell_face_idx[i_cell]; idx_face < pcell_face_idx[i_cell+1]; ++idx_face) {
        int i_face = PDM_ABS(pcell_face[idx_face])-1;
        for(int idx_vtx = pface_vtx_idx[i_face]; idx_vtx < pface_vtx_idx[i_face+1]; ++idx_vtx) {
          int i_vtx = PDM_ABS(pface_vtx[idx_vtx])-1;

          for (int i_dim = 0; i_dim < 3; i_dim++) {
            double x = vtx_coord[3*i_vtx + i_dim];

            if (x < _extents[i_dim]) {
              _extents[i_dim] = x;
            }
            if (x > _extents[3+i_dim]) {
              _extents[3+i_dim] = x;
            }
          }
        }
      }
    }

    free(pface_vtx);
  }


  const int dim = 3;
  double l_extents[6] = { HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                         -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};

  for(int i_part = 0; i_part < n_part; ++i_part) {
    for (int i = 0; i < pn_cell[i_part]; i++) {
      for (int k = 0; k < 3; k++) {
        l_extents[k    ] = PDM_MIN (l_extents[k    ], box_extents[i_part][6*i + k    ]);
        l_extents[k + 3] = PDM_MAX (l_extents[k + 3], box_extents[i_part][6*i + k + 3]);
      }
    }
  }

  double g_extents[6];
  PDM_MPI_Allreduce (l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce (l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX (max_range, g_extents[i+3] - g_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    g_extents[i]   -= max_range * 1.1e-3;
    g_extents[i+3] += max_range * 1.0e-3;
  }

  PDM_dbbtree_t* dbbt = PDM_dbbtree_create(comm, 3, g_extents);
  PDM_box_set_t *box_set = PDM_dbbtree_boxes_set (dbbt,
                                                  1,
                                                  pn_cell,
                                (const double **) box_extents,
                           (const PDM_g_num_t **) pcell_ln_to_gn);

  int          n_lines       = 0;
  double      *ray_coord     = NULL;
  PDM_g_num_t *pray_ln_to_gn = NULL;
  _create_wall_ray(comm,
                   n_part,
                   psurf_vtx,
                   psurf_face,
                   psurf_vtx_coord,
                   psurf_face_vtx_idx,
                   psurf_face_vtx,
                   psurf_face_ln_to_gn,
                   psurf_vtx_ln_to_gn,
                   &n_lines,
                   &pray_ln_to_gn,
                   &ray_coord);


  PDM_dbbtree_free (dbbt);
  PDM_box_set_destroy (&box_set);

  PDM_dist_cloud_surf_free(dist);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(psurf_vtx_coord    [i_part]);
    free(psurf_face_vtx_idx [i_part]);
    free(psurf_face_vtx     [i_part]);
    free(psurf_face_ln_to_gn[i_part]);
    free(psurf_vtx_ln_to_gn [i_part]);
    free(cell_center        [i_part]);
    free(velocity           [i_part]);
  }
  free(psurf_vtx );
  free(psurf_face);
  free(psurf_vtx_coord);
  free(psurf_face_vtx_idx);
  free(psurf_face_vtx);
  free(psurf_face_ln_to_gn);
  free(psurf_vtx_ln_to_gn);
  free(cell_center);
  free(pcell_ln_to_gn);
  free(pn_cell);
  free(velocity);

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
