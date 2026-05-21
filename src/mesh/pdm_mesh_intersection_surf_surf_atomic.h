#ifndef PDM_MESH_INTERSECTION_SURF_SURF_ATOMIC_H
#define PDM_MESH_INTERSECTION_SURF_SURF_ATOMIC_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"


/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct PDM_mesh_intersection_surf_surf_polygon_t {

  double *coord;
  int     n_edge;
  int    *face_vtx;
  int    *face_edge;
  int    *edge_vtx;

} PDM_mesh_intersection_surf_surf_polygon_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Compute the area and center of mass of the intersection between two polygons.
 *
 * \param [in]  poly_a    1st polygon
 * \param [in]  poly_b    2nd polygon
 * \param [out] area_ab   Surface area of intersection (signed)
 * \param [out] center_ab Center of mass of intersection (size = 3)
 */
void PDM_mesh_intersection_surf_surf_atomic_compute
(
  PDM_mesh_intersection_surf_surf_polygon_t *poly_a,
  PDM_mesh_intersection_surf_surf_polygon_t *poly_b,
  double                                    *area_ab,
  double                                    *center_ab
);

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_MESH_INTERSECTION_H */
