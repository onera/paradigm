/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2023       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_error.h"

#include "pdm_logging.h"

#include "pdm_vtk.h"
#include "pdm_reader_gamma.h"

#include "pdm_array.h"

#include "pdm_multipart.h"

#include "pdm_poly_vol_gen.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dcube_nodal_gen.h"

#include "pdm_isosurface_test_utils.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


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


static const int it_max = 40;
static double
_eval_mandelbrot
(
 const double  x,
 const double  y,
 const double  z
)
{
  PDM_UNUSED(z);
  int f = 0.;

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
  f = rk * log(PDM_MAX(1e-9, rk)) / PDM_MAX(1e-9, mag_dz);
  // if (it < it_max-1) {
  //   f = 1.;
  // } else {
  //   f = -1.;
  // }
  return f;
}


static
double
_eval_distance
(
  const double x,
  const double y,
  const double z,
  const double ctr[3]
)
{
  return 1.2*sqrt((x-ctr[0])*(x-ctr[0])+
                  (y-ctr[1])*(y-ctr[1])+
                  (z-ctr[2])*(z-ctr[2]));
}


static
double
_eval_cos
(
  const double x,
  const double y,
  const double z
)
{
  PDM_UNUSED(x);
  return cos(5.*(y + z));
}


static
void
_usage
(
  int exit_code
)
{
  printf
    ("\n"
     "  -h                             This message.\n\n"
     "  -n_part        <n>             Number of partitions (0 -> block-distributed).\n\n"
     "  -in            <filename>      Mesh file name (Gamma Mesh Format).\n\n"
     "  -sol           <filename>      Solution file name (Gamma Mesh Format).\n\n"
     "  -visu                          Enable exports for visualization.\n\n"
     "  -n_isovalues   <n>             Number of isovalues.\n\n"
     "  -isovalues     <v1 v2 ... vn>  Isovalues.\n\n"
     "  -elt_type      <t>             Volume element type (only for automatically generated mesh).\n\n"
     "  -randomize                     Enable mesh randomization (only for automatically generated mesh).\n\n"
     "  -n             <n>             Number of vertices per side (only for automatically generated mesh).\n\n"
     "  -use_part_mesh                 Use part_mesh structure (dmesh if n_part <= 0).\n\n"
     "  -edges                         Generate edges.\n\n"
     "  -local                         Build isosurface locally.\n\n"
     );

  exit(exit_code);
}

void
PDM_isosurface_test_utils_read_args
(
  int                    argc,
  char                 **argv,
  int                   *n_part,
  char                 **mesh_name,
  char                 **sol_name,
  int                   *visu,
  int                   *n_isovalues,
  double               **isovalues,
  PDM_Mesh_nodal_elt_t  *elt_type,
  int                   *randomize,
  PDM_g_num_t           *n_vtx_seg,
  int                   *use_part_mesh,
  int                   *generate_edges,
  int                   *local
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

    else if (strcmp(argv[i], "-in") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *mesh_name = argv[i];
      }
    }

    else if (strcmp(argv[i], "-sol") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *sol_name = argv[i];
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else if (strcmp(argv[i], "-n_isovalues") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_isovalues = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-isovalues") == 0) {
      PDM_malloc(*isovalues, *n_isovalues, double);
      for (int j = 0; j < *n_isovalues; j++) {
        i++;
        if (i >= argc)
          _usage(EXIT_FAILURE);
        else {
          (*isovalues)[j] = atof(argv[i]);
        }
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

    else if (strcmp(argv[i], "-rand") == 0) {
      *randomize = 1;
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

    else if (strcmp(argv[i], "-use_part_mesh") == 0) {
      *use_part_mesh = 1;
    }

    else if (strcmp(argv[i], "-edges") == 0) {
      *generate_edges = 1;
    }

    else if (strcmp(argv[i], "-local") == 0) {
      *local = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}


void
PDM_isosurface_test_utils_analytic_field_function
(
 const double  x,
 const double  y,
 const double  z,
 double       *value
)
{
  *value = cos(5*x)*cos(6*y)*cos(7*z);
}



void
PDM_isosurface_test_utils_gen_mesh
(
  PDM_MPI_Comm          comm,
  const char           *filename,
  int                   n_part,
  PDM_g_num_t           n_vtx_seg,
  int                   randomize,
  PDM_Mesh_nodal_elt_t  elt_type,
  int                   generate_edges,
  PDM_multipart_t     **mpart,
  PDM_part_mesh_t      *pmesh,
  PDM_dmesh_t         **out_dmesh
)
{
  int dim = PDM_Mesh_nodal_elt_dim_get(elt_type);

  if (n_part > 0) {
    int n_domain = 1;
    *mpart = PDM_multipart_create(n_domain,
                                  &n_part,
                                  PDM_FALSE,
                                  PDM_SPLIT_DUAL_WITH_HILBERT,
                                  PDM_PART_SIZE_HOMOGENEOUS,
                                  NULL,
                                  comm,
                                  PDM_OWNERSHIP_KEEP);
  }

  if (filename != NULL) {
    PDM_dmesh_nodal_t *dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                                          filename,
                                                          0,
                                                          0);

    if (0) {
      if (dim==3) {
        PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC,  "isosurface_3d_ngon_volume");
      }
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "isosurface_3d_ngon_surface");
    }

    if (n_part > 0) {
      PDM_multipart_dmesh_nodal_set(*mpart, 0, dmn);
      PDM_multipart_compute(*mpart);
    }
    else {
      PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1,
                                                                              comm,
                                                                              PDM_OWNERSHIP_USER);
      PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm,
                                               0,
                                               dmn);

      PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                       PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                       PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

      PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm,
                                         0,
                                         out_dmesh);

      PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);
      double *vtx_coord = PDM_DMesh_nodal_coord_get(dmn, PDM_OWNERSHIP_USER); // ugly but saving lives
      PDM_dmesh_vtx_coord_set(*out_dmesh, vtx_coord, PDM_OWNERSHIP_KEEP);
    }

    PDM_DMesh_nodal_free(dmn);
  }

  else if (elt_type == PDM_MESH_NODAL_POLY_3D) {
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
                     randomize,
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
    PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
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

    if (n_part > 0) {
      PDM_multipart_dmesh_set(*mpart, 0, dmesh);

      PDM_multipart_compute(*mpart);

      PDM_dmesh_free(dmesh);

      // Make edges great again (ugly AF)
      if (generate_edges) {
        int          *pn_face        = NULL;
        int          *pn_vtx         = NULL;
        int         **pface_vtx_idx  = NULL;
        int         **pface_vtx      = NULL;
        PDM_g_num_t **pface_ln_to_gn = NULL;
        PDM_g_num_t **pvtx_ln_to_gn  = NULL;
        PDM_malloc(pn_face       , n_part, int          );
        PDM_malloc(pn_vtx        , n_part, int          );
        PDM_malloc(pface_vtx_idx , n_part, int         *);
        PDM_malloc(pface_vtx     , n_part, int         *);
        PDM_malloc(pface_ln_to_gn, n_part, PDM_g_num_t *);
        PDM_malloc(pvtx_ln_to_gn , n_part, PDM_g_num_t *);

        int         **pface_edge_idx = NULL;
        int         **pface_edge     = NULL;
        int          *pn_edge        = NULL;
        int         **pedge_vtx      = NULL;
        PDM_g_num_t **pedge_ln_to_gn = NULL;

        for (int i_part = 0; i_part < n_part; i_part++) {
          PDM_multipart_part_connectivity_get(*mpart,
                                              0,
                                              i_part,
                                              PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                              &pface_vtx_idx[i_part],
                                              &pface_vtx    [i_part],
                                              PDM_OWNERSHIP_KEEP);

          pn_face[i_part] = PDM_multipart_part_ln_to_gn_get(*mpart,
                                                            0,
                                                            i_part,
                                                            PDM_MESH_ENTITY_FACE,
                                                            &pface_ln_to_gn[i_part],
                                                            PDM_OWNERSHIP_KEEP);

          pn_vtx[i_part] = PDM_multipart_part_ln_to_gn_get(*mpart,
                                                           0,
                                                           i_part,
                                                           PDM_MESH_ENTITY_VTX,
                                                           &pvtx_ln_to_gn[i_part],
                                                           PDM_OWNERSHIP_KEEP);
        }

        PDM_compute_face_edge_from_face_vtx(comm,
                                            n_part,
                                            pn_face,
                                            pn_vtx,
                                            pface_vtx_idx,
                                            pface_vtx,
                                            pface_ln_to_gn,
                                            pvtx_ln_to_gn,
                                            &pface_edge_idx,
                                            &pface_edge,
                                            &pn_edge,
                                            &pedge_vtx,
                                            &pedge_ln_to_gn);

        PDM_part_mesh_t *pm = (*mpart)->pmeshes[0].pmesh;
        for (int i_part = 0; i_part < n_part; i_part++) {
          PDM_part_mesh_connectivity_set(pm,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                         pface_edge    [i_part],
                                         pface_edge_idx[i_part],
                                         PDM_OWNERSHIP_KEEP);

          PDM_part_mesh_connectivity_set(pm,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                         pedge_vtx[i_part],
                                         NULL,
                                         PDM_OWNERSHIP_KEEP);

          PDM_part_mesh_n_entity_set(pm, i_part, PDM_MESH_ENTITY_EDGE, pn_edge[i_part]);
          PDM_part_mesh_entity_ln_to_gn_set(pm,
                                            i_part,
                                            PDM_MESH_ENTITY_EDGE,
                                            pedge_ln_to_gn[i_part],
                                            PDM_OWNERSHIP_KEEP);
        }
        PDM_free(pn_face       );
        PDM_free(pn_vtx        );
        PDM_free(pface_vtx_idx );
        PDM_free(pface_vtx     );
        PDM_free(pface_ln_to_gn);
        PDM_free(pvtx_ln_to_gn );
        PDM_free(pface_edge_idx);
        PDM_free(pface_edge    );
        PDM_free(pn_edge       );
        PDM_free(pedge_vtx     );
        PDM_free(pedge_ln_to_gn);
      }
    }
    else {
      *out_dmesh = dmesh;
    }
  }
  else {
    int n_vtx_seg_k = 1;
    if (dim==3) {
      n_vtx_seg_k = n_vtx_seg;
    }
    PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                          n_vtx_seg,
                                                          n_vtx_seg,
                                                          n_vtx_seg_k,
                                                          1,
                                                          0,
                                                          0,
                                                          0,
                                                          elt_type,
                                                          1,
                                                          PDM_OWNERSHIP_USER);

    PDM_dcube_nodal_gen_random_factor_set(dcube, (double) randomize);

    PDM_dcube_nodal_gen_build(dcube);

    PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

    PDM_dmesh_nodal_generate_distribution(dmn);

    if (n_part > 0) {
      PDM_multipart_dmesh_nodal_set(*mpart, 0, dmn);
      PDM_multipart_compute(*mpart);
    }
    else {
      PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1,
                                                                              comm,
                                                                              PDM_OWNERSHIP_USER);
      PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm,
                                               0,
                                               dmn);

      if (dim==3) {
        PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                         PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                         PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);
      } 
      else if (dim==2) {
        PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                         PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                         PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);        
      }

      PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm,
                                         0,
                                         out_dmesh);

      PDM_dmesh_compute_distributions(*out_dmesh);

      PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);
      double *vtx_coord = PDM_DMesh_nodal_coord_get(dmn, PDM_OWNERSHIP_USER); // ugly but saving lives
      PDM_dmesh_vtx_coord_set(*out_dmesh, vtx_coord, PDM_OWNERSHIP_KEEP);
    }
    PDM_DMesh_nodal_free(dmn);
    PDM_dcube_nodal_gen_free(dcube);
  }



  if (pmesh != NULL) {
    assert(n_part > 0);
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i_entity = PDM_MESH_ENTITY_CELL; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
        PDM_g_num_t *entity_ln_to_gn = NULL;
        int n_entity = PDM_multipart_part_ln_to_gn_get(*mpart,
                                                       0,
                                                       i_part,
                                                       i_entity,
                                                       &entity_ln_to_gn,
                                                       PDM_OWNERSHIP_KEEP);

        PDM_part_mesh_n_entity_set(pmesh, i_part, i_entity, n_entity);
        PDM_part_mesh_entity_ln_to_gn_set(pmesh,
                                          i_part,
                                          i_entity,
                                          entity_ln_to_gn,
                                          PDM_OWNERSHIP_USER);
      }

      int *cell_face_idx = NULL;
      int *cell_face     = NULL;
      PDM_multipart_part_connectivity_get(*mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                          &cell_face_idx,
                                          &cell_face,
                                          PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_connectivity_set(pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                     cell_face,
                                     cell_face_idx,
                                     PDM_OWNERSHIP_USER);

      if (generate_edges) {
        int *face_edge_idx = NULL;
        int *face_edge     = NULL;
        PDM_multipart_part_connectivity_get(*mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                            &face_edge_idx,
                                            &face_edge,
                                            PDM_OWNERSHIP_KEEP);
        assert(face_edge != NULL);
        PDM_part_mesh_connectivity_set(pmesh,
                                       i_part,
                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                       face_edge,
                                       face_edge_idx,
                                       PDM_OWNERSHIP_USER);

        int *edge_vtx_idx = NULL;
        int *edge_vtx     = NULL;
        PDM_multipart_part_connectivity_get(*mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                            &edge_vtx_idx,
                                            &edge_vtx,
                                            PDM_OWNERSHIP_KEEP);
        assert(edge_vtx != NULL);
        PDM_part_mesh_connectivity_set(pmesh,
                                       i_part,
                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                       edge_vtx,
                                       edge_vtx_idx,
                                       PDM_OWNERSHIP_USER);
      }
      else {
        if (dim==3) {
          int *face_vtx_idx = NULL;
          int *face_vtx     = NULL;
          PDM_multipart_part_connectivity_get(*mpart,
                                              0,
                                              i_part,
                                              PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                              &face_vtx_idx,
                                              &face_vtx,
                                              PDM_OWNERSHIP_KEEP);
          assert(face_vtx != NULL);
          PDM_part_mesh_connectivity_set(pmesh,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                         face_vtx,
                                         face_vtx_idx,
                                         PDM_OWNERSHIP_USER);
        }
      }

      double *vtx_coord = NULL;
      PDM_multipart_part_vtx_coord_get(*mpart, 0, i_part, &vtx_coord, PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_vtx_coord_set(pmesh, i_part, vtx_coord, PDM_OWNERSHIP_USER);


      int          n_surface             = 0;
      int         *surface_face_idx      = NULL;
      int         *surface_face          = NULL;
      PDM_g_num_t *surface_face_ln_to_gn = NULL;
      PDM_multipart_group_get(*mpart,
                              0,
                              i_part,
                              PDM_MESH_ENTITY_FACE,
                              &n_surface,
                              &surface_face_idx,
                              &surface_face,
                              &surface_face_ln_to_gn,
                              PDM_OWNERSHIP_KEEP);

      PDM_part_mesh_bound_concat_set(pmesh,
                                     i_part,
                                     PDM_BOUND_TYPE_FACE,
                                     n_surface,
                                     surface_face_idx,
                                     surface_face,
                                     surface_face_ln_to_gn,
                                     PDM_OWNERSHIP_USER);
    }
  }
}


void
PDM_isosurface_test_utils_gen_mesh_nodal
(
  PDM_MPI_Comm            comm,
  const char             *filename,
  int                     n_part,
  PDM_g_num_t             n_vtx_seg,
  int                     randomize,
  PDM_Mesh_nodal_elt_t    elt_type,
  PDM_part_mesh_nodal_t **out_pmn,
  PDM_dmesh_nodal_t     **out_dmn
)
{
  int dim = PDM_Mesh_nodal_elt_dim_get(elt_type);

  PDM_dmesh_nodal_t *dmn = NULL;

  if (filename != NULL) {
    dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                       filename,
                                       0,
                                       0);
  }

  else if (elt_type == PDM_MESH_NODAL_POLY_3D) {
    PDM_error(__FILE__, __LINE__, 0, "Poly3d not implemented yet\n");
  }

  else {
    int n_vtx_seg_k = 1;
    if (dim==3) {
      n_vtx_seg_k = n_vtx_seg;
    }
    PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                          n_vtx_seg,
                                                          n_vtx_seg,
                                                          n_vtx_seg_k,
                                                          1,
                                                          0,
                                                          0,
                                                          0,
                                                          elt_type,
                                                          1,
                                                          PDM_OWNERSHIP_USER);

    PDM_dcube_nodal_gen_random_factor_set(dcube, (double) randomize);

    PDM_dcube_nodal_gen_build(dcube);

    dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

    PDM_dcube_nodal_gen_free(dcube);
  }

  assert(dmn != NULL);

  PDM_dmesh_nodal_generate_distribution(dmn);

  if (n_part > 0) {
    int n_domain = 1;
    PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                  &n_part,
                                                  PDM_FALSE,
                                                  PDM_SPLIT_DUAL_WITH_HILBERT,
                                                  PDM_PART_SIZE_HOMOGENEOUS,
                                                  NULL,
                                                  comm,
                                                  PDM_OWNERSHIP_KEEP);

    PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
    PDM_multipart_compute(mpart);

    PDM_multipart_get_part_mesh_nodal(mpart,
                                      0,
                                      out_pmn,
                                      PDM_OWNERSHIP_USER);

    PDM_DMesh_nodal_free(dmn);
    PDM_multipart_free(mpart);
  }
  else {
    *out_dmn = dmn;
  }
}


void
PDM_isosurface_test_utils_compute_iso_field
(
  int     n_vtx,
  double *vtx_coord,
  double *vtx_field
)
{ 
  double ctr[3] = {0.1, 0.2, .3};
  double c2d[4] = {0.285, 0.01, 0., 0.};
  double c4d[4] = {-0.08, 0.0, -0.8, -0.03};
  for (int i = 0; i < 4; i++) { c4d[i] += 0.2; }

  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    double x = vtx_coord[3*i_vtx  ];
    double y = vtx_coord[3*i_vtx+1];
    double z = vtx_coord[3*i_vtx+2];
    int    ite = 0;
    if (1) {
      vtx_field[i_vtx] = _eval_distance(x,y,z,ctr)-0.3;
    }
    else {
      vtx_field[i_vtx] = _eval_julia4d(x,y,z,c2d,ite);
      vtx_field[i_vtx] = _eval_julia4d(x,y,z,c2d,ite);
      vtx_field[i_vtx] = _eval_julia4d(x,y,z,c4d,ite);
      vtx_field[i_vtx] = _eval_cos(x,y,z);
      vtx_field[i_vtx] = _eval_mandelbrot(x,y,z);
    }
  }
}

void
PDM_isosurface_test_utils_compute_itp_field
(
  int     n_vtx,
  double *vtx_coord,
  double *vtx_field
)
{
  double ctr[3] = {0.1, 0.2, .3};
  double c2d[4] = {0.285, 0.01, 0., 0.};
  double c4d[4] = {-0.08, 0.0, -0.8, -0.03};
  for (int i = 0; i < 4; i++) { c4d[i] += 0.2; }

  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    double x = vtx_coord[3*i_vtx  ];
    double y = vtx_coord[3*i_vtx+1];
    double z = vtx_coord[3*i_vtx+2];
    int    ite = 0;

    if (1) {
      vtx_field[i_vtx] = _eval_cos(x,y,z);
    }
    else {
      vtx_field[i_vtx] = _eval_julia4d(x,y,z,c2d,ite);
      vtx_field[i_vtx] = _eval_julia4d(x,y,z,c4d,ite);
      vtx_field[i_vtx] = _eval_distance(x,y,z,ctr);
      vtx_field[i_vtx] = _eval_mandelbrot(x,y,z);
    }
  }
}


void
PDM_isosurface_test_utils_dist_interpolation
(
  PDM_isosurface_t *isos,
  int               id_iso,
  double           *itp_dfield_vtx,
  double           *itp_dfield_face,
  double           *itp_dfield_cell,
  double          **iso_itp_dfield_vtx,
  double          **iso_itp_dfield_edge,
  double          **iso_itp_dfield_face
)
{
  int dim = isos->entry_mesh_dim;

  double *_iso_itp_dfield_vtx  = NULL;
  double *_iso_itp_dfield_edge = NULL;
  double *_iso_itp_dfield_face = NULL;

  /**
   * Vertex interpolation
   */
  PDM_part_to_part_t *ptp_vtx = NULL;
  PDM_isosurface_part_to_part_get(isos, id_iso, PDM_MESH_ENTITY_VTX,
                                 &ptp_vtx, PDM_OWNERSHIP_BAD_VALUE); 

  // > Exchange isofield
  double **recv_vtx_field = NULL;
  int request_vtx = -1;
  PDM_part_to_part_reverse_iexch(ptp_vtx,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(double),
                                 NULL,
                (const void  **)&itp_dfield_vtx,
                                 NULL,
                (      void ***)&recv_vtx_field,
                                &request_vtx);

  PDM_part_to_part_reverse_iexch_wait(ptp_vtx, request_vtx);

  int  *n_ref_lnum2 = NULL;
  int **  ref_lnum2 = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp_vtx, &n_ref_lnum2, &ref_lnum2);
  PDM_log_trace_array_int(ref_lnum2[0], n_ref_lnum2[0], "ref_lnum2");

  int         *vtx_dparent_idx  = NULL;
  double      *vtx_dparent_wght = NULL;
  int iso_dn_vtx = PDM_isosurface_dparent_weight_get(isos, id_iso, PDM_MESH_ENTITY_VTX,
                                                    &vtx_dparent_idx, &vtx_dparent_wght,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_calloc(_iso_itp_dfield_vtx, iso_dn_vtx, double);
  for (int i_vtx = 0; i_vtx < iso_dn_vtx; i_vtx++) {
    for (int i = vtx_dparent_idx[i_vtx]; i < vtx_dparent_idx[i_vtx+1]; i++) {
      _iso_itp_dfield_vtx[i_vtx] += vtx_dparent_wght[i]*recv_vtx_field[0][i];
    }
  }

  PDM_free(recv_vtx_field[0]);
  PDM_free(recv_vtx_field);

  *iso_itp_dfield_vtx  = _iso_itp_dfield_vtx;


  /**
   * Edge interpolation
   */
  PDM_part_to_part_t *ptp_edge = NULL;
  PDM_isosurface_part_to_part_get(isos, id_iso, PDM_MESH_ENTITY_EDGE,
                                 &ptp_edge, PDM_OWNERSHIP_BAD_VALUE); 

  // > Exchange isofield
  double **recv_edge_field = NULL;
  int request_edge = -1;
  PDM_part_to_part_reverse_iexch(ptp_edge,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(double),
                                 NULL,
                (const void  **)&itp_dfield_face,
                                 NULL,
                (      void ***)&recv_edge_field,
                                &request_edge);

  PDM_part_to_part_reverse_iexch_wait(ptp_edge, request_edge);

  n_ref_lnum2 = NULL;
    ref_lnum2 = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp_edge, &n_ref_lnum2, &ref_lnum2);
  PDM_log_trace_array_int(ref_lnum2[0], n_ref_lnum2[0], "ref_lnum2");

  int         *edge_dparent_idx  = NULL;
  double      *edge_dparent_wght = NULL;
  int iso_dn_edge = PDM_isosurface_dparent_weight_get(isos, id_iso, PDM_MESH_ENTITY_EDGE,
                                                    &edge_dparent_idx, &edge_dparent_wght,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_calloc(_iso_itp_dfield_edge, iso_dn_edge, double);
  for (int i_edge = 0; i_edge < iso_dn_edge; i_edge++) {
    for (int i = edge_dparent_idx[i_edge]; i < edge_dparent_idx[i_edge+1]; i++) {
      _iso_itp_dfield_edge[i_edge] += edge_dparent_wght[i]*recv_edge_field[0][i];
    }
  }

  PDM_free(recv_edge_field[0]);
  PDM_free(recv_edge_field);

  *iso_itp_dfield_edge = _iso_itp_dfield_edge;


  /**
   * Face interpolation
   */
  if (dim==3) {
    PDM_part_to_part_t *ptp_face = NULL;
    PDM_isosurface_part_to_part_get(isos, id_iso, PDM_MESH_ENTITY_FACE,
                                   &ptp_face, PDM_OWNERSHIP_BAD_VALUE); 

    // > Exchange isofield
    double **recv_face_field = NULL;
    int request_face = -1;
    PDM_part_to_part_reverse_iexch(ptp_face,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                   1,
                                   sizeof(double),
                                   NULL,
                  (const void  **)&itp_dfield_cell,
                                   NULL,
                  (      void ***)&recv_face_field,
                                  &request_face);

    PDM_part_to_part_reverse_iexch_wait(ptp_face, request_face);

    n_ref_lnum2 = NULL;
      ref_lnum2 = NULL;
    PDM_part_to_part_ref_lnum2_get(ptp_face, &n_ref_lnum2, &ref_lnum2);
    PDM_log_trace_array_int(ref_lnum2[0], n_ref_lnum2[0], "ref_lnum2");

    int         *face_dparent_idx  = NULL;
    double      *face_dparent_wght = NULL;
    int iso_dn_face = PDM_isosurface_dparent_weight_get(isos, id_iso, PDM_MESH_ENTITY_FACE,
                                                      &face_dparent_idx, &face_dparent_wght,
                                                       PDM_OWNERSHIP_KEEP);

    PDM_calloc(_iso_itp_dfield_face, iso_dn_face, double);
    for (int i_face = 0; i_face < iso_dn_face; i_face++) {
      for (int i = face_dparent_idx[i_face]; i < face_dparent_idx[i_face+1]; i++) {
        _iso_itp_dfield_face[i_face] += face_dparent_wght[i]*recv_face_field[0][i];
      }
    }

    PDM_free(recv_face_field[0]);
    PDM_free(recv_face_field);
    
    *iso_itp_dfield_face = _iso_itp_dfield_face;
  }
}


void
PDM_isosurface_test_utils_part_interpolation
(
  PDM_isosurface_t *isos,
  int               id_iso,
  int               n_part,
  int               local,
  double          **itp_field_vtx,
  double          **itp_field_face,
  double          **itp_field_cell,
  double         ***iso_itp_field_vtx,
  double         ***iso_itp_field_edge,
  double         ***iso_itp_field_face
)
{
  int dim = isos->entry_mesh_dim;

  // TODO: interpolate face and cell fields
  PDM_UNUSED(itp_field_face);
  PDM_UNUSED(itp_field_cell);
  PDM_UNUSED(iso_itp_field_edge);
  PDM_UNUSED(iso_itp_field_face);

  double **_iso_itp_field_vtx  = NULL;
  double **_iso_itp_field_edge = NULL;
  double **_iso_itp_field_face = NULL;

  PDM_malloc(_iso_itp_field_vtx , n_part, double *);
  PDM_malloc(_iso_itp_field_edge, n_part, double *);
  if (dim==3) {
    PDM_malloc(_iso_itp_field_face, n_part, double *);
  }

  if (local) {
    // Local    
    for (int i_part = 0; i_part < n_part; i_part++) {
      int    *iso_vtx_parent_idx;
      int    *iso_vtx_parent;
      double *iso_vtx_parent_wght;
      int iso_n_vtx = PDM_isosurface_parent_weight_get(isos,
                                                       id_iso,
                                                       i_part,
                                                       PDM_MESH_ENTITY_VTX,
                                                       &iso_vtx_parent_idx,
                                                       &iso_vtx_parent_wght,
                                                       PDM_OWNERSHIP_KEEP);

      PDM_isosurface_local_parent_get(isos,
                                      id_iso,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &iso_vtx_parent_idx,
                                      &iso_vtx_parent,
                                      PDM_OWNERSHIP_KEEP);

      PDM_malloc(_iso_itp_field_vtx[i_part], iso_n_vtx, double);
      for (int i_vtx = 0; i_vtx < iso_n_vtx; i_vtx++) {
        _iso_itp_field_vtx[i_part][i_vtx] = 0.;
        for (int i = iso_vtx_parent_idx[i_vtx]; i < iso_vtx_parent_idx[i_vtx+1]; i++) {
          int i_parent = iso_vtx_parent[i] - 1;
          _iso_itp_field_vtx[i_part][i_vtx] += iso_vtx_parent_wght[i] * itp_field_vtx[i_part][i_parent];
        }
      }

      int    *iso_edge_parent_idx;
      int    *iso_edge_parent;
      double *iso_edge_parent_wght;
      int iso_n_edge = PDM_isosurface_local_parent_get(isos,
                                                       id_iso,
                                                       i_part,
                                                       PDM_MESH_ENTITY_EDGE,
                                                       &iso_edge_parent_idx,
                                                       &iso_edge_parent,
                                                       PDM_OWNERSHIP_KEEP);

      iso_n_edge = PDM_isosurface_parent_weight_get(isos,
                                                    id_iso,
                                                    i_part,
                                                    PDM_MESH_ENTITY_EDGE,
                                                    &iso_edge_parent_idx,
                                                    &iso_edge_parent_wght,
                                                    PDM_OWNERSHIP_KEEP);
      PDM_malloc(_iso_itp_field_edge[i_part], iso_n_edge, double);
      for (int i_edge = 0; i_edge < iso_n_edge; i_edge++) {
        _iso_itp_field_edge[i_part][i_edge] = 0.;
        for (int i = iso_edge_parent_idx[i_edge]; i < iso_edge_parent_idx[i_edge+1]; i++) {
          int i_parent = iso_edge_parent[i] - 1;
          _iso_itp_field_edge[i_part][i_edge] += iso_edge_parent_wght[i] * itp_field_face[i_part][i_parent];
        }
      }

      if (dim==3) {
        int    *iso_face_parent_idx;
        int    *iso_face_parent;
        double *iso_face_parent_wght;
        int iso_n_face = PDM_isosurface_local_parent_get(isos,
                                                         id_iso,
                                                         i_part,
                                                         PDM_MESH_ENTITY_FACE,
                                                         &iso_face_parent_idx,
                                                         &iso_face_parent,
                                                         PDM_OWNERSHIP_KEEP);

        iso_n_face = PDM_isosurface_parent_weight_get(isos,
                                                      id_iso,
                                                      i_part,
                                                      PDM_MESH_ENTITY_FACE,
                                                      &iso_face_parent_idx,
                                                      &iso_face_parent_wght,
                                                      PDM_OWNERSHIP_KEEP);

        PDM_malloc(_iso_itp_field_face[i_part], iso_n_face, double);
        for (int i_face = 0; i_face < iso_n_face; i_face++) {
          _iso_itp_field_face[i_part][i_face] = 0.;
          for (int i = iso_face_parent_idx[i_face]; i < iso_face_parent_idx[i_face+1]; i++) {
            int i_parent = iso_face_parent[i] - 1;
            _iso_itp_field_face[i_part][i_face] += iso_face_parent_wght[i] * itp_field_cell[i_part][i_parent];
          }
        }
      }
    }

  } // End if LOCAL

  else {
    // Reequilibrate

    // Interpolate a vtx-based field
    PDM_part_to_part_t *ptp_vtx = NULL;
    PDM_isosurface_part_to_part_get(isos,
                                    id_iso,
                                    PDM_MESH_ENTITY_VTX,
                                    &ptp_vtx,
                                    PDM_OWNERSHIP_KEEP);

    double **recv_vtx_field = NULL;
    int request_vtx = -1;
    PDM_part_to_part_reverse_iexch(ptp_vtx,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                   1,
                                   sizeof(double),
                                   NULL,
                  (const void  **) itp_field_vtx,
                                   NULL,
                  (      void ***) &recv_vtx_field,
                                   &request_vtx);

    PDM_part_to_part_reverse_iexch_wait(ptp_vtx, request_vtx);

    for (int i_part = 0; i_part < n_part; i_part++) {
      int    *iso_vtx_parent_idx;
      double *iso_vtx_parent_wght;
      int iso_n_vtx = PDM_isosurface_parent_weight_get(isos, id_iso, i_part, PDM_MESH_ENTITY_VTX,
                                                      &iso_vtx_parent_idx,
                                                      &iso_vtx_parent_wght,
                                                       PDM_OWNERSHIP_KEEP);

      PDM_malloc(_iso_itp_field_vtx[i_part], iso_n_vtx, double);
      for (int i_vtx = 0; i_vtx < iso_n_vtx; i_vtx++) {
        _iso_itp_field_vtx[i_part][i_vtx] = 0.;
        for (int i = iso_vtx_parent_idx[i_vtx]; i < iso_vtx_parent_idx[i_vtx+1]; i++) {
          _iso_itp_field_vtx[i_part][i_vtx] += iso_vtx_parent_wght[i] * recv_vtx_field[i_part][i];
        }
      }

      PDM_free(recv_vtx_field[i_part]);
    } // End loop on parts
    PDM_free(recv_vtx_field);



    // Exchange parent face gnums just for fun
    PDM_part_to_part_t *ptp_edge = NULL;
    PDM_isosurface_part_to_part_get(isos,
                                    id_iso,
                                    PDM_MESH_ENTITY_EDGE,
                                    &ptp_edge,
                                    PDM_OWNERSHIP_KEEP);

    double **recv_edge_field = NULL;
    int request_edge = -1;
    PDM_part_to_part_reverse_iexch(ptp_edge,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                   1,
                                   sizeof(double),
                                   NULL,
                  (const void  **) itp_field_face,
                                   NULL,
                  (      void ***) &recv_edge_field,
                                   &request_edge);

    PDM_part_to_part_reverse_iexch_wait(ptp_edge, request_edge);

    // convert gnums to doubles for vtk export
    for (int i_part = 0; i_part < n_part; i_part++) {
      int    *iso_edge_parent_idx;
      double *iso_edge_parent_wght;
      int iso_n_edge = PDM_isosurface_parent_weight_get(isos, id_iso, i_part, PDM_MESH_ENTITY_EDGE,
                                                      &iso_edge_parent_idx,
                                                      &iso_edge_parent_wght,
                                                       PDM_OWNERSHIP_KEEP);

      PDM_malloc(_iso_itp_field_edge[i_part], iso_n_edge, double);
      for (int i_edge = 0; i_edge < iso_n_edge; i_edge++) {
        _iso_itp_field_edge[i_part][i_edge] = 0.;
        for (int i = iso_edge_parent_idx[i_edge]; i < iso_edge_parent_idx[i_edge+1]; i++) {
          _iso_itp_field_edge[i_part][i_edge] += iso_edge_parent_wght[i] * recv_edge_field[i_part][i];
        }
      }

      PDM_free(recv_edge_field[i_part]);
    }
    PDM_free(recv_edge_field);


    // Exchange parent cell gnums just for fun
    if (dim==3) {
      PDM_part_to_part_t *ptp_face = NULL;
      PDM_isosurface_part_to_part_get(isos,
                                      id_iso,
                                      PDM_MESH_ENTITY_FACE,
                                      &ptp_face,
                                      PDM_OWNERSHIP_KEEP);


      double **recv_face_field = NULL;
      int request_face = -1;
      PDM_part_to_part_reverse_iexch(ptp_face,
                                     PDM_MPI_COMM_KIND_P2P,
                                     PDM_STRIDE_CST_INTERLACED,
                                     PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                     1,
                                     sizeof(double),
                                     NULL,
                    (const void  **) itp_field_cell,
                                     NULL,
                    (      void ***) &recv_face_field,
                                     &request_face);

      PDM_part_to_part_reverse_iexch_wait(ptp_face, request_face);

      // convert gnums to doubles for vtk export
      for (int i_part = 0; i_part < n_part; i_part++) {
        int    *iso_face_parent_idx;
        double *iso_face_parent_wght;
        int iso_n_face = PDM_isosurface_parent_weight_get(isos, id_iso, i_part, PDM_MESH_ENTITY_FACE,
                                                        &iso_face_parent_idx,
                                                        &iso_face_parent_wght,
                                                         PDM_OWNERSHIP_KEEP);

        PDM_malloc(_iso_itp_field_face[i_part], iso_n_face, double);
        for (int i_face = 0; i_face < iso_n_face; i_face++) {
          _iso_itp_field_face[i_part][i_face] = 0.;
          for (int i = iso_face_parent_idx[i_face]; i < iso_face_parent_idx[i_face+1]; i++) {
            _iso_itp_field_face[i_part][i_face] += iso_face_parent_wght[i] * recv_face_field[i_part][i];
          }
        }

        PDM_free(recv_face_field[i_part]);
      }
      PDM_free(recv_face_field);
    }

  } // End if REEQUILIBRATE

  *iso_itp_field_vtx  = _iso_itp_field_vtx;
  *iso_itp_field_edge = _iso_itp_field_edge;
  if (dim==3) {
    *iso_itp_field_face = _iso_itp_field_face;
  }
}


void
PDM_isosurface_test_utils_dist_vtk
(
  PDM_isosurface_t *isos,
  int               id_iso,
  const double     *iso_vtx_fld,
  const double     *iso_edge_fld,
  const double     *iso_face_fld,
  PDM_MPI_Comm      comm
)
{
  int  debug = 1;
  char out_name[999];

  int dim = isos->entry_mesh_dim;

  // > Vertices
  double *dvtx_coords = NULL;
  int     dn_vtx = PDM_isosurface_dvtx_coord_get(isos, id_iso, &dvtx_coords, PDM_OWNERSHIP_KEEP);
  if (1) {
    log_trace("dn_vtx = %d\n", dn_vtx);
    PDM_log_trace_array_double(dvtx_coords, 3*dn_vtx, "dvtx_coords ::");
  }

  // > Edges
  int         *dedge_vtx_idx = NULL;
  PDM_g_num_t *dedge_vtx     = NULL;
  PDM_g_num_t dn_edge = PDM_isosurface_dconnectivity_get(isos, id_iso, 
                                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                        &dedge_vtx_idx,
                                                        &dedge_vtx,
                                                         PDM_OWNERSHIP_KEEP);
  if (debug) {
    log_trace("dn_edge = "PDM_FMT_G_NUM"\n", dn_edge);
    PDM_log_trace_array_long(dedge_vtx, 2*dn_edge, "dedge_vtx ::");
  }

  // > Faces
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t dn_face        = 0;
  if (dim==3) {
    dn_face = PDM_isosurface_dconnectivity_get(isos, id_iso, 
                                               PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                              &dface_vtx_idx,
                                              &dface_vtx,
                                               PDM_OWNERSHIP_KEEP);
    if (debug) {
      log_trace("dn_face = "PDM_FMT_G_NUM"\n", dn_face);
      int size_face_vtx = dface_vtx_idx[dn_face];
      PDM_log_trace_array_int (dface_vtx_idx, dn_face      , "dface_vtx_idx ::");
      PDM_log_trace_array_long(dface_vtx    , size_face_vtx, "dface_vtx     ::");
    }
  }

  // > Groups
  int          n_group         = 0;
  int         *dedge_group_idx = NULL;
  PDM_g_num_t *dedge_group     = NULL;
  if (dim==3) {
    n_group = PDM_isosurface_dgroup_get(isos, id_iso, 
                                        PDM_MESH_ENTITY_EDGE,
                                       &dedge_group_idx,
                                       &dedge_group,
                                        PDM_OWNERSHIP_KEEP);
    if (debug) {
      log_trace("n_group = %d\n", n_group);
      int size_edge_group = dedge_group_idx[n_group];
      PDM_log_trace_array_int (dedge_group_idx, n_group        , "dedge_group_idx ::");
      PDM_log_trace_array_long(dedge_group    , size_edge_group, "dedge_group ::");
    }
  }


  /**
   * Prepare dmesh_nodal for ouput
   */
  PDM_g_num_t  n_entity[3] = {dn_vtx, dn_edge, dn_face};
  PDM_g_num_t gn_entity[3];
  PDM_MPI_Allreduce(n_entity, gn_entity, 3, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_log_trace_array_long(gn_entity, 3, "gn_entity ::");
  PDM_dmesh_nodal_t *iso_dmn = PDM_DMesh_nodal_create(comm,
                                                      dim,
                                                      gn_entity[0],
                                                      0,
                                                      gn_entity[2],
                                                      gn_entity[1]);

  PDM_DMesh_nodal_coord_set(iso_dmn, dn_vtx, dvtx_coords, PDM_OWNERSHIP_USER);

  int edge_section = PDM_DMesh_nodal_section_add(iso_dmn,
                                                 PDM_GEOMETRY_KIND_RIDGE,
                                                 PDM_MESH_NODAL_BAR2);
  PDM_DMesh_nodal_section_std_set(iso_dmn,
                                  PDM_GEOMETRY_KIND_RIDGE,
                                  edge_section,
                                  dn_edge,
                                  dedge_vtx,
                                  PDM_OWNERSHIP_USER);

  if (dim==3) {
    int face_section = PDM_DMesh_nodal_section_add(iso_dmn,
                                                   PDM_GEOMETRY_KIND_SURFACIC,
                                                   PDM_MESH_NODAL_POLY_2D);

    PDM_DMesh_nodal_section_poly2d_set(iso_dmn,
                                       PDM_GEOMETRY_KIND_SURFACIC,
                                       face_section,
                                       dn_face,
                                       dface_vtx_idx,
                                       dface_vtx,
                                       PDM_OWNERSHIP_USER);

    PDM_DMesh_nodal_section_group_elmt_set(iso_dmn,
                                           PDM_GEOMETRY_KIND_RIDGE,
                                           n_group,
                                           dedge_group_idx,
                                           dedge_group,
                                           PDM_OWNERSHIP_USER);
  }

  sprintf(out_name, "iso_edge_id_%d.vtk", id_iso);
  PDM_dmesh_nodal_dump_vtk(iso_dmn, PDM_GEOMETRY_KIND_RIDGE, out_name);

  const char *fld_name_vtx[]  = {"itp_vtx_fld"};
  const char *fld_name_edge[] = {"itp_edge_fld"};
  PDM_dmesh_nodal_dump_vtk_with_field(iso_dmn, PDM_GEOMETRY_KIND_RIDGE, out_name,
                                      1,
                                      fld_name_vtx,
                                      &iso_vtx_fld,
                                      1,
                                      fld_name_edge,
                                      &iso_edge_fld);
  if (dim==3) {
    const char *fld_name_face[] = {"itp_face_fld"};
    sprintf(out_name, "iso_face_id_%d.vtk", id_iso);
    PDM_dmesh_nodal_dump_vtk_with_field(iso_dmn, PDM_GEOMETRY_KIND_SURFACIC, out_name,
                                        1,
                                        fld_name_vtx,
                                        &iso_vtx_fld,
                                        1,
                                        fld_name_face,
                                        &iso_face_fld);
  }

  PDM_DMesh_nodal_free(iso_dmn);
}


void
PDM_isosurface_test_utils_part_vtk
(
  PDM_isosurface_t   *isos,
  int                 id_iso,
  int                 n_part,
  double           **iso_vtx_fld,
  double           **iso_edge_fld,
  double           **iso_face_fld,
  PDM_MPI_Comm        comm
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int dim = isos->entry_mesh_dim;
  
  for (int i_part=0; i_part<n_part; ++i_part) {

    /**
     * Get isosurface data
     */
    // > Vertex
    double      *iso_vtx_coord = NULL;
    PDM_g_num_t *iso_vtx_gnum  = NULL;
    int iso_n_vtx = PDM_isosurface_vtx_coord_get(isos, id_iso, i_part, &iso_vtx_coord, PDM_OWNERSHIP_KEEP);
    PDM_isosurface_ln_to_gn_get(isos, id_iso, i_part, PDM_MESH_ENTITY_VTX, &iso_vtx_gnum, PDM_OWNERSHIP_KEEP);

    // > Compute isovalue tag
    int    *vtx_isovalue_idx = NULL;
    double *vtx_isovalue_tag = NULL;
    int _n_isovalues = PDM_isosurface_isovalue_entity_idx_get(isos, id_iso, i_part, PDM_MESH_ENTITY_VTX,
                                                             &vtx_isovalue_idx, PDM_OWNERSHIP_KEEP);
    PDM_malloc(vtx_isovalue_tag, iso_n_vtx, double);
    for (int i_isovalue=0; i_isovalue<_n_isovalues; ++i_isovalue) {
      for (int i_entity=vtx_isovalue_idx[i_isovalue  ];
               i_entity<vtx_isovalue_idx[i_isovalue+1]; ++i_entity) {
        vtx_isovalue_tag[i_entity] = (double) (i_isovalue+1);
      }
    }


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

    int          n_group             = 0;
    int         *group_edge_idx      = NULL;
    int         *group_edge          = NULL;
    PDM_g_num_t *group_edge_ln_to_gn = NULL;
    double      *fld_group_edge      = NULL;
    PDM_malloc(fld_group_edge, iso_n_edge, double);
    n_group = PDM_isosurface_group_get(isos,
                                       id_iso,
                                       i_part,
                                       PDM_MESH_ENTITY_EDGE,
                                       &group_edge_idx,
                                       &group_edge,
                                       &group_edge_ln_to_gn,
                                       PDM_OWNERSHIP_KEEP);
    for (int i_group = 0; i_group < n_group; i_group++) {
      for (int i = group_edge_idx[i_group]; i < group_edge_idx[i_group+1]; i++) {
        fld_group_edge[group_edge[i]-1] = (double) (i_group + 1);
      }
    }

    // > Convert group into tag
    double *iso_edge_tag1      = NULL;
    double *iso_edge_tag2      = NULL;
    PDM_calloc(iso_edge_tag1     , iso_n_edge, double);
    PDM_calloc(iso_edge_tag2     , iso_n_edge, double);
    if (dim==3) {
      iso_edge_n_group = PDM_isosurface_group_get(isos, id_iso, i_part, PDM_MESH_ENTITY_EDGE,
                                                 &iso_edge_group_idx, &iso_edge_group_lnum, &iso_edge_group_gnum,
                                                  PDM_OWNERSHIP_KEEP);
      for (int i_group=0; i_group<iso_edge_n_group; ++i_group) {
        int i_beg_group = iso_edge_group_idx[i_group  ];
        int i_end_group = iso_edge_group_idx[i_group+1];
        for (int i_read=i_beg_group; i_read<i_end_group; ++i_read) {
          int         edge_lnum = iso_edge_group_lnum[i_read];
          if (iso_edge_tag1[edge_lnum-1]<1.e-9) {
            iso_edge_tag1     [edge_lnum-1] = (double) (i_group+1);
          }
          iso_edge_tag2     [edge_lnum-1] = (double) (i_group+1);
        }
      }
    }

    // > Compute isovalue tag
    int    *edge_isovalue_idx = NULL;
    double *edge_isovalue_tag = NULL;
    _n_isovalues = PDM_isosurface_isovalue_entity_idx_get(isos, id_iso, i_part, PDM_MESH_ENTITY_EDGE,
                                                         &edge_isovalue_idx, PDM_OWNERSHIP_KEEP);
    PDM_malloc(edge_isovalue_tag, iso_n_edge, double);
    for (int i_isovalue=0; i_isovalue<_n_isovalues; ++i_isovalue) {
      for (int i_entity=edge_isovalue_idx[i_isovalue  ];
               i_entity<edge_isovalue_idx[i_isovalue+1]; ++i_entity) {
        edge_isovalue_tag[i_entity] = (double) (i_isovalue+1);
      }
    }


    // > Face
    int          iso_n_face        = 0;
    int         *iso_face_vtx_idx  = NULL;
    int         *iso_face_vtx      = NULL;
    PDM_g_num_t *iso_face_gnum     = NULL;
    double      *face_isovalue_tag = NULL;
    if (dim==3) {
      iso_n_face = PDM_isosurface_connectivity_get(isos, id_iso, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX, &iso_face_vtx_idx, &iso_face_vtx, PDM_OWNERSHIP_KEEP);
      PDM_isosurface_ln_to_gn_get(isos, id_iso, i_part, PDM_MESH_ENTITY_FACE, &iso_face_gnum, PDM_OWNERSHIP_KEEP);

      int *isovalue_face_idx = NULL;
      _n_isovalues = PDM_isosurface_isovalue_entity_idx_get(isos,
                                                            id_iso,
                                                            i_part,
                                                            PDM_MESH_ENTITY_FACE,
                                                            &isovalue_face_idx,
                                                            PDM_OWNERSHIP_KEEP);

      PDM_malloc(face_isovalue_tag, iso_n_face, double);
      for (int i_isovalue = 0; i_isovalue < _n_isovalues; i_isovalue++) {
        for (int i_face = isovalue_face_idx[i_isovalue]; i_face < isovalue_face_idx[i_isovalue+1]; i_face++) {
          face_isovalue_tag[i_face] = i_isovalue+1;
        }
      }
    }

    /**
     * Write isosurface
     */
    char filename[999];


    // > 3D ngon old version
    sprintf(filename, "iso_edge_id_%i_part_%i_rank_%i.vtk", id_iso, i_part, i_rank);

    const char   *vtx_field_name [] = {"i_isovalue", "itp_field"};
    const double *vtx_field_value[] = {vtx_isovalue_tag, (const double *) iso_vtx_fld[i_part]};

    const char   *edge_field_name [] = {"i_isovalue", "i_group1", "i_group2", "itp_edge_fld"};
    const double *edge_field_value[] = {edge_isovalue_tag, iso_edge_tag1, iso_edge_tag2, iso_edge_fld[i_part]};

    PDM_vtk_write_std_elements_ho_with_vtx_field(filename,
                                                 1,
                                                 iso_n_vtx,
                                                 iso_vtx_coord,
                                                 iso_vtx_gnum,
                                                 PDM_MESH_NODAL_BAR2,
                                                 iso_n_edge,
                                                 iso_edge_vtx,
                                                 iso_edge_gnum,
                                                 4,
                                                 edge_field_name,
                                                 edge_field_value,
                                                 2,
                                                 vtx_field_name,
                                                 vtx_field_value);

    if (dim==3) {
      sprintf(filename, "iso_face_id_%i_part_%i_rank_%i.vtk", id_iso, i_part, i_rank);

      const char   *face_field_name [] = {"i_isovalue", "itp_field"};
      const double *face_field_value[] = {face_isovalue_tag, iso_face_fld[i_part]};

      int *elt_type = PDM_array_const_int(iso_n_face, PDM_MESH_NODAL_POLY_2D);
      PDM_vtk_write_unstructured_grid(filename,
                                      iso_n_vtx,
                                      iso_vtx_coord,
                                      iso_vtx_gnum,
                                      iso_n_face,
            (PDM_Mesh_nodal_elt_t *)  elt_type,
                                      iso_face_vtx_idx,
                                      iso_face_vtx,
                                      iso_face_gnum,
                                      2,
                                      face_field_name,
                                      face_field_value,
                                      2,
                                      vtx_field_name,
                                      vtx_field_value);
      PDM_free(elt_type);

    }
    PDM_free( vtx_isovalue_tag);
    PDM_free(edge_isovalue_tag);
    PDM_free(face_isovalue_tag);

    PDM_free(iso_edge_tag1);
    PDM_free(iso_edge_tag2);
    PDM_free(fld_group_edge);
  }

}


#ifdef  __cplusplus
}
#endif
