/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_mesh_intersection_surf_surf_atomic.h"
#include "pdm.h"
#include "pdm_error.h"
#include "pdm_line.h"
#include "pdm_logging.h"
#include "pdm_priv.h"

/*============================================================================
 * Global variables
 *============================================================================*/

// static int dbg_enabled = 0;

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static inline void
_get_edge_vtx
(
  PDM_mesh_intersection_surf_surf_polygon_t  *poly,
  int                                         i,
  int                                        *i_vtx0,
  int                                        *i_vtx1,
  double                                    **p0,
  double                                    **p1
)
{
  if (poly->face_vtx != NULL) {
    *i_vtx0 = poly->face_vtx[ i                ] - 1;
    *i_vtx1 = poly->face_vtx[(i+1)%poly->n_edge] - 1;
  }
  else {
    int i_edge = poly->face_edge[i];
    if (i_edge < 0) {
      i_edge = -i_edge - 1;
      *i_vtx0 = poly->edge_vtx[2*i_edge+1] - 1;
      *i_vtx1 = poly->edge_vtx[2*i_edge  ] - 1;
    }
    else {
      i_edge = i_edge - 1;
      *i_vtx0 = poly->edge_vtx[2*i_edge  ] - 1;
      *i_vtx1 = poly->edge_vtx[2*i_edge+1] - 1;
    }
  }

  *p0 = &poly->coord[3*(*i_vtx0)];
  *p1 = &poly->coord[3*(*i_vtx1)];
}


static int
_clip_segment
(
  double *uv0,
  double *uv1,
  double *uv2
)
{
  double f0, f1, t, s;

  int n_comp = 0;

  // Clip segment to keep intersection with the positive quadrant (u >= 0, v >= 0)
  for (int i = 0; i < 2; i++) {

    f0 = uv0[i];
    f1 = uv1[i];

    if (f0*f1 < 0) {
      t = f0/(f0 - f1);

      if (f0 < 0) {
        // keep segment [t, 1]
        uv0[0] = (1-t)*uv0[0] + t*uv1[0];
        uv0[1] = (1-t)*uv0[1] + t*uv1[1];
      }
      else { // f1 < 0
        // keep segment [0, t]
        uv1[0] = (1-t)*uv0[0] + t*uv1[0];
        uv1[1] = (1-t)*uv0[1] + t*uv1[1];
      }
    }
    else if (f0 < 0) { // f1 < 0 as well
      // The segment does not intersect the positive quadrant
      return n_comp;
    }
  }

  // Find intersection with diagonal (u + v = 1)
  f0 = 1 - uv0[0] - uv0[1];
  f1 = 1 - uv1[0] - uv1[1];
  if (f0*f1 < 0) {
    n_comp = 2;
    t = f0/(f0 - f1);
    uv2[0] = (1-t)*uv0[0] + t*uv1[0];
    uv2[1] = (1-t)*uv0[1] + t*uv1[1];
  }
  else {
    // no intersection
    n_comp = 1;
  }

  // Project outer part on diagonal
  if (f0 < 0) {
    s = 1./(1 - f0);
    uv0[0] *= s;
    uv0[1] *= s;
  }
  if (f1 < 0) {
    s = 1./(1 - f1);
    uv1[0] *= s;
    uv1[1] *= s;
  }

  return n_comp;
}


#define ONE_THIRD 0.3333333333333333
static inline void
_geom_comp
(
  const double *uv0,
  const double *uv1,
        double *area,
        double *center
)
{
  *area = uv0[0]*uv1[1] - uv0[1]*uv1[0];
  center[0] = ONE_THIRD*(uv0[0] + uv1[0]);
  center[1] = ONE_THIRD*(uv0[1] + uv1[1]);
}
#undef ONE_THIRD


static inline void
_vector_ab
(
        double ab[3],
  const double a[3],
  const double b[3]
)
{
  ab[0] = b[0] - a[0];
  ab[1] = b[1] - a[1];
  ab[2] = b[2] - a[2];
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

void PDM_mesh_intersection_surf_surf_atomic_compute
(
  PDM_mesh_intersection_surf_surf_polygon_t *poly_a,
  PDM_mesh_intersection_surf_surf_polygon_t *poly_b,
  double                                    *area_ab,
  double                                    *center_ab
)
{
  *area_ab = 0.;
  for (int i = 0; i < 3; i++) {
    center_ab[i] = 0.;
  }

  int i_vtx_ref = -1;

  double *pi, *pj;

  double ai[3] = {0., 0., 0.};
  double aj[3], bi[3], bj[3];
  double uvi[2], uvj[2], uvk[2];

  for (int ia = 0; ia < poly_a->n_edge; ia++) {

    int i_vtx_a, j_vtx_a;
    _get_edge_vtx(poly_a, ia, &i_vtx_a, &j_vtx_a, &pi, &pj);

    i_vtx_ref = (i_vtx_ref < 0) ? i_vtx_a : i_vtx_ref;

    if ((i_vtx_a == i_vtx_ref) || (j_vtx_a == i_vtx_ref)) {
      continue;
    }

    _vector_ab(ai, &poly_a->coord[3*i_vtx_ref], pi);
    _vector_ab(aj, &poly_a->coord[3*i_vtx_ref], pj);

    double aiai = PDM_DOT_PRODUCT(ai, ai);
    double aiaj = PDM_DOT_PRODUCT(ai, aj);
    double ajaj = PDM_DOT_PRODUCT(aj, aj);

    double det_aij = aiai*ajaj - aiaj*aiaj;

    if (det_aij <= 0) {
      // vertices i_vtx_ref, i_vtx_a and j_vtx_a are collinear => ignore this subtriangle
      continue;
    }

    double inv_det_aij = 1./det_aij;

    double normal_aij[3];
    PDM_CROSS_PRODUCT(normal_aij, ai, aj);

    double area_aij = 0.5*PDM_MODULE(normal_aij);

    for (int ib = 0; ib < poly_b->n_edge; ib++) {

      int i_vtx_b, j_vtx_b;
      _get_edge_vtx(poly_b, ib, &i_vtx_b, &j_vtx_b, &pi, &pj);

      _vector_ab(bi, &poly_a->coord[3*i_vtx_ref], pi);
      _vector_ab(bj, &poly_a->coord[3*i_vtx_ref], pj);


      double biai = PDM_DOT_PRODUCT(bi, ai);
      double biaj = PDM_DOT_PRODUCT(bi, aj);
      double bjai = PDM_DOT_PRODUCT(bj, ai);
      double bjaj = PDM_DOT_PRODUCT(bj, aj);

      uvi[0] = biai*ajaj - biaj*aiaj;
      uvj[0] = bjai*ajaj - bjaj*aiaj;

      if ((uvi[0] <= 0) && (uvj[0] <= 0)) {
        continue;
      }

      uvi[1] = aiai*biaj - aiaj*biai;
      uvj[1] = aiai*bjaj - aiaj*bjai;

      if ((uvi[1] <= 0) && (uvj[1] <= 0)) {
        continue;
      }

      uvi[0] *= inv_det_aij;
      uvi[1] *= inv_det_aij;
      uvj[0] *= inv_det_aij;
      uvj[1] *= inv_det_aij;


      int n_comp = _clip_segment(uvi, uvj, uvk);

      if (n_comp == 0) {
        continue;
      }

      double area_uv;
      double center_uv[2] = {0., 0.};
      if (n_comp == 1) {
        _geom_comp(uvi, uvj, &area_uv, center_uv);
      }
      else { // n_comp == 2
        double area_uv_ik, area_uv_kj;
        double center_uv_ik[2], center_uv_kj[2];
        _geom_comp(uvi, uvk, &area_uv_ik, center_uv_ik);
        _geom_comp(uvk, uvj, &area_uv_kj, center_uv_kj);
        area_uv = area_uv_ik + area_uv_kj;
PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
        if (area_uv != 0) {
          double inv_area_uv = 1./area_uv;
          for (int l = 0; l < 2; l++) {
            center_uv[l] = (area_uv_ik*center_uv_ik[l] + area_uv_kj*center_uv_kj[l]) * inv_area_uv;
          }
        }
PDM_GCC_SUPPRESS_WARNING_POP
      }

      area_uv *= area_aij;

      *area_ab += area_uv;
      for (int l = 0; l < 3; l++) {
        center_ab[l] += area_uv * (
          poly_a->coord[3*i_vtx_ref+l] +
          center_uv[0]*ai[l] +
          center_uv[1]*aj[l]
        );
      }

    } // End loop on edges B
  } // End loop on edges A

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (*area_ab != 0) {
    double inv_area_ab = 1./(*area_ab);
    for (int l = 0; l < 3; l++) {
      center_ab[l] *= inv_area_ab;
    }
  }
PDM_GCC_SUPPRESS_WARNING_POP
}


