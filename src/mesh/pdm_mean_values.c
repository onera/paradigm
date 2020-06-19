/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_plane.h"
#include "pdm_line.h"
#include "pdm_polygon.h"
#include "pdm_geom_elem.h"
#include "pdm_triangulate.h"

#include "fvmc_defs.h"
#include "fvmc_triangulate.h"
#include "fvmc_point_location.h"
#include "fvmc_ho_basis.h"

#include "pdm_mean_values.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Private function definition
 *============================================================================*/

static void
computePolygonMeanValues(const int           n_dist_points,
                         const fvmc_lnum_t  *dist_locations,
                         const fvmc_coord_t *dist_coords,
                         const int          *meshConnectivityIndex,
                         const int          *meshConnectivity,
                         const double       *meshVertexCoords,
                         /* const std::vector <int>& nDistBarCoords, */
                         /* std::vector <double>& distBarCoords */
                         const int nDistBarCoords[],
                         double    *distBarCoords)

{
  /* Boucle sur les points distants */

  fvmc_coord_t coo_point_dist[3];

  /* Tableaux locaux */

  const double eps_base = 1e-10;
  double *coo_som_fac = NULL;
  double *s           = NULL;
  double *dist        = NULL;
  double *aire        = NULL;
  double *proScal     = NULL;
  int size = 0;

  for (int ipoint =  0; ipoint < n_dist_points; ipoint++ ) {

    double *_distBarCoords = distBarCoords + nDistBarCoords[ipoint];

    /* Initialisation - Copie locale */

    int isOnEdge = 0;
    int isVertex = 0;
    int ielt = dist_locations[ipoint] - 1;

    int nbr_som_fac =  meshConnectivityIndex[ielt+1] -
                       meshConnectivityIndex[ielt];
    coo_point_dist[0] = dist_coords[3*ipoint];
    coo_point_dist[1] = dist_coords[3*ipoint + 1];
    coo_point_dist[2] = dist_coords[3*ipoint + 2];

    if (ipoint == 0) {
      size = nbr_som_fac;

      coo_som_fac = malloc (sizeof(double) * size * 3);
      s           = malloc (sizeof(double) * size * 3);
      dist        = malloc (sizeof(double) * size);
      aire        = malloc (sizeof(double) * size);
      proScal     = malloc (sizeof(double) * size);
    }
    else {
      if (size < nbr_som_fac) {
        size = nbr_som_fac;

        coo_som_fac = realloc (coo_som_fac, sizeof(double) * size * 3);
        s           = realloc (s,           sizeof(double) * size * 3);
        dist        = realloc (dist,        sizeof(double) * size);
        aire        = realloc (aire,        sizeof(double) * size);
        proScal     = realloc (proScal,     sizeof(double) * size);
      }
    }

    for (int isom = 0; isom < nbr_som_fac; isom++) {
      coo_som_fac[3*isom]   =
        meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)];

      coo_som_fac[3*isom+1] =
        meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)+1];

      coo_som_fac[3*isom+2] =
        meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)+2];
    }

    /* Projection sur un plan moyen */

    double bary[3];
    PDM_polygon_compute_barycenter (nbr_som_fac, &(coo_som_fac[0]), bary);

    double n[3] = {0, 0, 1};
    PDM_plane_normal (nbr_som_fac, &(coo_som_fac[0]), n);

    PDM_plane_projection2 (coo_point_dist, bary, n, coo_point_dist);

    for (int isom = 0; isom < nbr_som_fac; isom++) {

      double *pt1 = &(coo_som_fac[0]) + 3 *isom;
      PDM_plane_projection2 (pt1, bary, n, pt1);

    }

    double bounds[6] = {DBL_MAX, -DBL_MAX,
                        DBL_MAX, -DBL_MAX,
                        DBL_MAX, -DBL_MAX};

    for (int isom = 0; isom < nbr_som_fac; isom++) {
      bounds[0] = PDM_MIN (bounds[0], coo_som_fac[3*isom]);
      bounds[1] = PDM_MAX (bounds[1], coo_som_fac[3*isom]);

      bounds[2] = PDM_MIN (bounds[2], coo_som_fac[3*isom + 1]);
      bounds[3] = PDM_MAX (bounds[3], coo_som_fac[3*isom + 1]);

      bounds[4] = PDM_MIN (bounds[4], coo_som_fac[3*isom + 2]);
      bounds[5] = PDM_MAX (bounds[5], coo_som_fac[3*isom + 2]);
    }


    /* Verification que le point est dans l'element */

    if (PDM_polygon_point_in (coo_point_dist,
                              nbr_som_fac,
                              &(coo_som_fac[0]),
                              bounds,
                              n) != 1) {

      double closestPoint[3];
      double dist_min = DBL_MAX;

      for (int k = 0; k < nbr_som_fac; k++) {
        double *p1 = &(coo_som_fac[3 * k]);
        double *p2 = &(coo_som_fac[3 * ((k+1) % nbr_som_fac)]);
        double closest[3];
        double t;

        double dist2 = PDM_line_distance (coo_point_dist,
                                          p1,
                                          p2,
                                          &t,
                                          closest);

        if (dist2 < dist_min) {
          dist_min = dist2;
          closestPoint[0] = closest[0];
          closestPoint[1] = closest[1];
          closestPoint[2] = closest[2];
        }
      }

      coo_point_dist[0] = closestPoint[0];
      coo_point_dist[1] = closestPoint[1];
      coo_point_dist[2] = closestPoint[2];

    }

    /* Calcul des coordonnnees barycentriques */

    double min_dist = DBL_MAX;
    for (int isom = 0; isom < nbr_som_fac; isom++) {

      int inext = (isom + 1) % nbr_som_fac;
      double *vect = &s[0] + 3*isom;
      double l_edge;
      vect[0] = coo_som_fac[3*inext]   - coo_som_fac[3*isom];
      vect[1] = coo_som_fac[3*inext+1] - coo_som_fac[3*isom+1];
      vect[2] = coo_som_fac[3*inext+2] - coo_som_fac[3*isom+2];
      l_edge  = PDM_MODULE (vect);

      min_dist = PDM_MIN (l_edge, min_dist);
    }

    double eps = PDM_MAX (min_dist * eps_base, 1.e-30);

    for (int isom = 0; isom < nbr_som_fac; isom++) {

      double *vect = &s[0] + 3*isom;
      vect[0] = coo_som_fac[3*isom]   - coo_point_dist[0];
      vect[1] = coo_som_fac[3*isom+1] - coo_point_dist[1];
      vect[2] = coo_som_fac[3*isom+2] - coo_point_dist[2];
      dist[isom] = PDM_MODULE (vect);

    }

    int currentVertex;
    for (int isom = 0; isom < nbr_som_fac; isom++) {
      int inext = (isom + 1) % nbr_som_fac;
      double *vect1 = &s[0] + 3 * isom;
      double *vect2 = &s[0] + 3 * inext;
      double pvect[3];

      proScal[isom] = PDM_DOT_PRODUCT (vect1, vect2);
      PDM_CROSS_PRODUCT(pvect, vect1, vect2);

      double sign = PDM_DOT_PRODUCT (pvect, n);
      aire[isom] = PDM_MODULE(pvect);

      if (sign < 0) {
        aire[isom] = -aire[isom];
      }

      if (dist[isom] <= eps) {

        isVertex = 1;
        currentVertex = isom;
        break;
      }

      else if ((fabs(aire[isom]) <= eps)  && (proScal[isom] < 0)) {

        isOnEdge = 1;
        currentVertex = isom;
        break;

      }

    }

    /* Le point distant est un sommet */

    if (isVertex) {
      for (int isom = 0; isom < nbr_som_fac; isom++)
        _distBarCoords[isom] = 0.;
      _distBarCoords[currentVertex] = 1.;
    }

    /* Le point distant est sur arete */

    else if (isOnEdge) {

      for (int isom = 0; isom < nbr_som_fac; isom++)
        _distBarCoords[isom] = 0.;

      int nextPoint = (currentVertex + 1) % nbr_som_fac;

      _distBarCoords[currentVertex] =
        dist[nextPoint]     / (dist[nextPoint]+dist[currentVertex]);
      _distBarCoords[nextPoint]     =
        dist[currentVertex] / (dist[nextPoint]+dist[currentVertex]);

    }

    /* Cas general */

    else {

      double sigma = 0;
      for (int isom = 0; isom < nbr_som_fac; isom++) {
        double coef = 0.;
        int previousVertex = (isom - 1 + nbr_som_fac) % nbr_som_fac;
        int nextVertex = (isom + 1) % nbr_som_fac;

        if (fabs(aire[previousVertex]) > eps)
          coef += (dist[previousVertex] - proScal[previousVertex]/dist[isom]) / aire[previousVertex];
        if (fabs(aire[isom]) > eps)
          coef += (dist[nextVertex]     - proScal[isom]/dist[isom])           / aire[isom];
        sigma += coef;
        _distBarCoords[isom] = coef;

      }

      if (fabs(sigma) >= eps ) {
        for (int isom = 0; isom < nbr_som_fac; isom++) {
          _distBarCoords[isom] /= sigma;
        }
      }

      else {

        double abs_sigma = fabs(sigma);
        printf("Warning : Mise à NAN %f %f\n", abs_sigma,  eps);
        for (int isom = 0; isom < nbr_som_fac; isom++) {
          _distBarCoords[isom] = NAN;
        }
      }

      /* Check Result */

      for (int isom = 0; isom <  nbr_som_fac; isom++) {
        if ( _distBarCoords[isom] != _distBarCoords[isom] ||
             _distBarCoords[isom] < 0. ||
             _distBarCoords[isom] > 1. ) {

          double dist_min = DBL_MAX;
          int k_min = 0;
          double t_min;

          for (int k = 0; k < nbr_som_fac; k++) {
            _distBarCoords[k] = 0.0;
          }

          for (int k = 0; k < nbr_som_fac; k++) {
            double *p1 = &(coo_som_fac[3 * k]);
            double *p2 = &(coo_som_fac[3 * ((k+1) % nbr_som_fac)]);
            double closest[3];
            double t;

            double dist2 = PDM_line_distance (coo_point_dist,
                                              p1,
                                              p2,
                                              &t,
                                              closest);
            if (dist2 < dist_min) {
              t_min = t;
              k_min = k;
            }
          }

          _distBarCoords[k_min] = 1 - t_min;
          _distBarCoords[(k_min + 1) % nbr_som_fac] = t_min;

          break;

        }

      }

    }

    if (0 == 1) {
      if ((n_dist_points == 1) && (dist_locations[0] == 1)) {

        PDM_printf("coord %i %i :", ipoint+1, ielt+1);
        PDM_printf(" %12.5e %12.5e %12.5e", dist_coords[3*ipoint],
                    dist_coords[3*ipoint+1],
                    dist_coords[3*ipoint+2] );
        PDM_printf("\n");

        PDM_printf("coo b %i :", ipoint+1);
        for (int isom = 0; isom < nbr_som_fac; isom++) {
          PDM_printf(" %f", _distBarCoords[isom]);
        }
        PDM_printf("\n");
      }
    }
  }

  free (coo_som_fac);
  free (s);
  free (aire);
  free (dist);
  free (proScal);
}


//-->> duplicates from pdm_geom_elem.c --> change scope?
static const double GEOM_EPS_VOL  = 1e-9; /*!< Constant value used to compute geomtric epsilon for volume */

static const double GEOM_EPS_SURF = 1e-9; /*!< Constant value used to compute geomtric epsilon for surface */
//<<--

static void
compute3DMeanValuesPoly(const double point_coords[],
                        const int    n_poly_faces,
                        const int    n_poly_vertex,
                        const int    faceDirection[],
                        const int    faceToVertexIdx[],
                        const int    faceToVertex[],
                        const double vertex_coords[],
                        const double characteristicLength,
                        const float  distElt,
                        double       distBarCoords[])
{

  //
  // Polyhedron
  //


  // Mise a jour des tableaux locaux

  /* double *coo_som_face = malloc (sizeof(double) * 3 * n_poly_vertex); */
  double *dist = malloc (sizeof(double) * n_poly_vertex);
  double *s = malloc (sizeof(double) * 3 * n_poly_vertex);
  double angle[3]; //angle[  v(i) v v(i+1) ]
  double normale[9]; //normale

  double sigma = 0;

  int isOnFace = 0;

  double eps_loc = PDM_geom_elem_geometric_epsilon (characteristicLength, GEOM_EPS_VOL);

  /**** Inialisation du tableau des coordonnees temporaires a 0 ****/

  for (int isom = 0; isom < n_poly_vertex; isom++)
    distBarCoords[isom] = 0.;

  for (int isom = 0; isom < n_poly_vertex; isom++) {
    s[3 * isom    ] = vertex_coords[3 * isom    ] - point_coords[0];
    s[3 * isom + 1] = vertex_coords[3 * isom + 1] - point_coords[1];
    s[3 * isom + 2] = vertex_coords[3 * isom + 2] - point_coords[2];
  }

  if (distElt > 1.) {

    //
    // Search clostest face
    //

    int n_vtx_max = 0;
    for (int i = 0; i < n_poly_faces; i++) {
      const int n_vtx = faceToVertexIdx[i+1] - faceToVertexIdx[i];
      if (n_vtx > n_vtx_max) {
        n_vtx_max = n_vtx;
      }
    }

    double * poly_vertex = (double *) malloc (3 * sizeof(double) * n_vtx_max);

    double dist_min = DBL_MAX;
    int iface = -1;

    for (int i = 0; i < n_poly_faces; i++) {
      const int n_vtx = faceToVertexIdx[i+1] - faceToVertexIdx[i];
      int k = 0;

      for (int j = faceToVertexIdx[i]; j < faceToVertexIdx[i+1]; j++) {
        int i_vertex = faceToVertex[j] - 1;
        poly_vertex[k++] = vertex_coords[3 * i_vertex];
        poly_vertex[k++] = vertex_coords[3 * i_vertex + 1];
        poly_vertex[k++] = vertex_coords[3 * i_vertex + 2];

      }

      double closest[3];
      double dist_face;

      PDM_polygon_evaluate_position (point_coords,
                                     n_vtx, poly_vertex, closest,
                                     &dist_face);

      if (dist_face < dist_min) {
        dist_min = dist_face;
        iface = i;
      }

    }

    free (poly_vertex);

    //
    // Compute mean value for this face if point is on face
    //

    if (iface >= 0) {

      const int n_face_vertex = faceToVertexIdx[iface + 1]
                              - faceToVertexIdx[iface];

      //
      // Copy

      int distBarCoordsFaceIdx[2];
      distBarCoordsFaceIdx[0] = 0;
      distBarCoordsFaceIdx[1] = n_face_vertex;
      double *distBarCoordsFace = malloc (sizeof(double) * n_face_vertex);
      int face_location = iface + 1;

      /*computePolygonMeanValues(1,
                               &face_location,
                               point_coords,
                               faceToVertexIdx,
                               faceToVertex,
                               vertex_coords,
                               distBarCoordsFaceIdx,
                               distBarCoordsFace);*/
      PDM_geom_elem_compute_polygon_barycentric_coordinates (1,
                                                             &face_location,
                                                             point_coords,
                                                             faceToVertexIdx,
                                                             faceToVertex,
                                                             vertex_coords,
                                                             distBarCoordsFaceIdx,
                                                             distBarCoordsFace);

      for (int j = 0; j < n_poly_vertex; j++)
        distBarCoords[j] = 0.;

      for (int j = 0; j < n_face_vertex; j++) {
        int vertex = faceToVertex[faceToVertexIdx[iface]+j] - 1;
        distBarCoords[vertex] = distBarCoordsFace[j];
      }

      /* distBarCoordsFaceIdx.clear(); */
      free (distBarCoordsFace);

    }

    else {

      //
      // If point is not in a face, get the closest vertex
      //

      double normS = sqrt(s[3*0]*s[3*0] + s[3*0+1]*s[3*0+1] + s[3*0+2]*s[3*0+2]);
      int closestVertex = 0;
      for (int isom = 1; isom < n_poly_vertex; isom++) {

        double nextNormS = sqrt (s[3*isom]*s[3*isom]
                                 + s[3*isom+1]*s[3*isom+1]
                                 + s[3*isom+2]*s[3*isom+2]);
        if (nextNormS < normS) {
          closestVertex = isom;
          normS = nextNormS;
        }
      }

      distBarCoords[closestVertex] = 1;

    }

  }

  else if (distElt > 0 ) {

    //
    // Check if point is on a face
    //

    const int   n_points  = 1;
    fvmc_lnum_t point_ids = 0;

    //
    // Search clostest face
    //

    fvmc_lnum_t face_location = -1;
    float face_distance = -1.;

    fvmc_point_dist_closest_polygon(3,
                                    n_poly_faces,
                                    faceToVertexIdx,
                                    faceToVertex,
                                    vertex_coords,
                                    n_points,
                                    &point_ids,
                                    point_coords,
                                    &face_location,
                                    &face_distance);

    //
    // Compute mean value for this face if point is on face
    //

    if (face_location > 0 && face_distance < eps_loc) {

      isOnFace = 1;

      const int face_location_idx = face_location - 1;
      const int n_face_vertex = faceToVertexIdx[face_location_idx+1]
                              - faceToVertexIdx[face_location_idx];

      //
      // Copy

      int distBarCoordsFaceIdx[2];
      distBarCoordsFaceIdx[0] = 0;
      distBarCoordsFaceIdx[1] = n_face_vertex;
      double *distBarCoordsFace = malloc (sizeof(double) * n_face_vertex);

      /*computePolygonMeanValues(1,
                               &face_location,
                               point_coords,
                               faceToVertexIdx,
                               faceToVertex,
                               vertex_coords,
                               distBarCoordsFaceIdx,
                               distBarCoordsFace);*/
      PDM_geom_elem_compute_polygon_barycentric_coordinates (1,
                                                             &face_location,
                                                             point_coords,
                                                             faceToVertexIdx,
                                                             faceToVertex,
                                                             vertex_coords,
                                                             distBarCoordsFaceIdx,
                                                             distBarCoordsFace);

      for (int j = 0; j < n_poly_vertex; j++)
        distBarCoords[j] = 0.;

      for (int j = 0; j < n_face_vertex; j++) {
        int vertex = faceToVertex[faceToVertexIdx[face_location_idx]+j] - 1;

        distBarCoords[vertex] = distBarCoordsFace[j];

        /* distBarCoordsFaceIdx.clear(); */
        free (distBarCoordsFace);

      }

    }

    //
    // General alogorithm for point in polyhedron
    //

    if (!isOnFace) {

      for (int isom = 0; isom < n_poly_vertex; isom++) {

        dist[isom] = sqrt(s[3*isom    ] * s[3*isom    ]
                          + s[3*isom + 1] * s[3*isom + 1]
                          + s[3*isom + 2] * s[3*isom + 2]);

        s[3*isom]     /= dist[isom];
        s[3*isom + 1] /= dist[isom];
        s[3*isom + 2] /= dist[isom];

      }

      //
      // Second loop on faces to commpute barycentric coordinates
      //


      int s_triangle_vertices = 9;
      int *triangle_vertices = malloc (sizeof(int) * s_triangle_vertices); //Nombre de sommets apres decoupage en triangle

      for (int iface = 0; iface < n_poly_faces; iface++) {

        const int n_vertex_fac = faceToVertexIdx[iface + 1]
                               - faceToVertexIdx[iface];

        const int ind_fac_som = faceToVertexIdx[iface];

        fvmc_triangulate_state_t *fvmc_triangulate = fvmc_triangulate_state_create(n_vertex_fac);

        int triangle_vertice_size = (n_vertex_fac-2) * 3;

        if (s_triangle_vertices < triangle_vertice_size) {
          s_triangle_vertices = triangle_vertice_size;
          triangle_vertices = realloc (triangle_vertices, sizeof(int) * s_triangle_vertices);
        }

        //
        // Face triangulation
        //

        int n_triangles;

        if (n_vertex_fac == 4) {

          n_triangles = fvmc_triangulate_quadrangle(3,
                                                    vertex_coords,
                                                    NULL,
                                                    faceToVertex + ind_fac_som,
                                                    &(triangle_vertices[0]));

        }

        else if (n_vertex_fac > 4) {

          n_triangles = fvmc_triangulate_polygon(3,
                                                 n_vertex_fac,
                                                 vertex_coords,
                                                 NULL,
                                                 faceToVertex + ind_fac_som,
                                                 FVMC_TRIANGULATE_MESH_DEF,
                                                 &(triangle_vertices[0]),
                                                 fvmc_triangulate);
        }

        else {
          n_triangles = 1;
          for (int i = 0; i < 3; i++) {
            triangle_vertices[i] = faceToVertex[ind_fac_som + i];
          }
        }

        //
        // Loop on triangles
        //

        for (int itri = 0; itri < n_triangles; itri++) {

          //
          // Check triangle surface
          //

          const int i = triangle_vertices[3*itri    ] - 1;
          int j, k;
          if (faceDirection[iface] < 0) {
            j = triangle_vertices[3*itri + 1] - 1;
            k = triangle_vertices[3*itri + 2] - 1;
          }
          else {
            j = triangle_vertices[3*itri + 2] - 1;
            k = triangle_vertices[3*itri + 1] - 1;
          }

          const double coo_ijx = vertex_coords[3*j]   - vertex_coords[3*i];
          const double coo_ijy = vertex_coords[3*j+1] - vertex_coords[3*i+1];
          const double coo_ijz = vertex_coords[3*j+2] - vertex_coords[3*i+2];
          const double coo_ikx = vertex_coords[3*k]   - vertex_coords[3*i];
          const double coo_iky = vertex_coords[3*k+1] - vertex_coords[3*i+1];
          const double coo_ikz = vertex_coords[3*k+2] - vertex_coords[3*i+2];

          const double areaTri_ijk = 0.5 * sqrt((coo_ijy * coo_ikz - coo_ijz * coo_iky)
                                                * (coo_ijy * coo_ikz - coo_ijz * coo_iky)
                                                + (coo_ijz * coo_ikx - coo_ijx * coo_ikz)
                                                * (coo_ijz * coo_ikx - coo_ijx * coo_ikz)
                                                + (coo_ijx * coo_iky - coo_ijy * coo_ikx)
                                                * (coo_ijx * coo_iky - coo_ijy * coo_ikx));

          double eps_face = PDM_geom_elem_geometric_epsilon (characteristicLength, GEOM_EPS_SURF);

          if (fabs(areaTri_ijk) > eps_face) {

            for (int isom = 0; isom < 3; isom++) {

              int isuiv;
              int iprec;
              double prod_scal;
              double mod;

              if (faceDirection[iface] < 0) {
                iprec = triangle_vertices[3*itri + (isom + 2) % 3] - 1;
                isuiv = triangle_vertices[3*itri + (isom + 1) % 3] - 1;
              }
              else {
                iprec = triangle_vertices[3*itri + (isom + 1) % 3] - 1;
                isuiv = triangle_vertices[3*itri + (isom + 2) % 3] - 1;
              }

              prod_scal = s[3*iprec    ] * s[3*isuiv    ]
                        + s[3*iprec + 1] * s[3*isuiv + 1]
                        + s[3*iprec + 2] * s[3*isuiv + 2];

              angle[isom] = acos(prod_scal); //s

              normale[3 * isom    ] =  s[3*iprec + 1] * s[3*isuiv + 2]
                                     - s[3*iprec + 2] * s[3*isuiv + 1];
              normale[3 * isom + 1] =  s[3*iprec + 2] * s[3*isuiv    ]
                                     - s[3*iprec    ] * s[3*isuiv + 2];
              normale[3 * isom + 2] =  s[3*iprec    ] * s[3*isuiv + 1]
                                     - s[3*iprec + 1] * s[3*isuiv    ];

              /// verifier norm

              mod = sqrt(normale[3*isom    ] * normale[3*isom    ]
                         + normale[3*isom + 1] * normale[3*isom + 1]
                         + normale[3*isom + 2] * normale[3*isom + 2]);

              if (mod <  eps_face) {
                normale[3*isom    ] = 0.;
                normale[3*isom + 1] = 0.;
                normale[3*isom + 2] = 0.;
              }

              else {

                normale[3*isom    ] /= mod;
                normale[3*isom + 1] /= mod;
                normale[3*isom + 2] /= mod;
              }
            }

            for (int isom = 0; isom < 3; isom++) {

              double ps_nij_njk; //a ameliorer
              double ps_nki_njk; //a ameliorer
              double ps_ei_njk;  //a ameliorer

              const int iprec = (isom + 2) % 3;
              const int isuiv = (isom + 1) % 3;

              ps_nij_njk = normale[3 * isom    ] * normale[3 * isuiv    ]
                + normale[3 * isom + 1] * normale[3 * isuiv + 1]
                + normale[3 * isom + 2] * normale[3 * isuiv + 2];

              ps_nki_njk = normale[3 * isom    ] * normale[3 * iprec    ]
                + normale[3 * isom + 1] * normale[3 * iprec + 1]
                + normale[3 * isom + 2] * normale[3 * iprec + 2];

              // ps_ei_njk --> sur la face


              const int ivertex_tri = triangle_vertices[3*itri + isom] - 1;
              ps_ei_njk = s[3*ivertex_tri    ] * normale[3*isom]
                        + s[3*ivertex_tri + 1] * normale[3*isom + 1]
                        + s[3*ivertex_tri + 2] * normale[3*isom + 2];

              // vérifier ps_ei_njk

              if (fabs(ps_ei_njk) >  eps_face) {
                distBarCoords[ivertex_tri] +=
                  (angle[isom] + angle[isuiv] * ps_nij_njk + angle[iprec] * ps_nki_njk)
                  / (2 * ps_ei_njk);
               }

            } // Loop on vertices

          } // Good triangle

        } // Loop on triangles

        fvmc_triangulate_state_destroy(fvmc_triangulate);

      } // Loop on faces

      for (int isom = 0; isom < n_poly_vertex; isom++) {

        distBarCoords[isom] /= dist[isom];
        sigma += distBarCoords[isom];

      }

      for (int isom = 0; isom < n_poly_vertex; isom++)
        distBarCoords[isom] = distBarCoords[isom] / sigma;

    } // End of general algorithm (if (!isonface))

    //
    // Output results
    //

    if (0 == 1) {

      double test[3];

      for (int i = 0; i < 3; i++)
        test[i] = 0;

      for (int isom = 0; isom < n_poly_vertex; isom++){

        test[0] += distBarCoords[isom] * vertex_coords[3*isom];
        test[1] += distBarCoords[isom] * vertex_coords[3*isom + 1];
        test[2] += distBarCoords[isom] * vertex_coords[3*isom + 2];

      }

      PDM_printf("point distant | verification \n");

      double dd = 0;
      for (int i = 0; i < 3; i++) {
        PDM_printf("  %f       |    %f \n",point_coords[i],test[i]);
        dd += (point_coords[i] - test[i]) * (point_coords[i] - test[i]);
      }

      if (sqrt(dd) > 1e-3)
        PDM_printf(" !!!! Erreur sur les coordonnees baryc directionf: %12.5e %i !!!!\n",sqrt(dd), isOnFace);
      else
        PDM_printf(" ++++ ok                                         : %12.5e %i ++++\n",sqrt(dd), isOnFace);

      PDM_printf("coord :");
      PDM_printf(" %12.5e %12.5e %12.5e", point_coords[0],
                  point_coords[1],
                  point_coords[2] );
      PDM_printf("\n");

      PDM_printf("coo b :");
      for (int isom = 0; isom < n_poly_vertex; isom++)
        PDM_printf(" %f", distBarCoords[isom]);

      PDM_printf("\n");

    }

  }

}



/*=============================================================================
 * Public function definition
 *============================================================================*/

/**
 * \brief Compute mean values of a list of points in a polygon
 *
 * \param [in]    n_pts            Number of points to locate
 * \param [in]    pts              xyz-Coordinates of points locate
 * \param [in]    n_vtx            Number of polygon vertices
 * \param [in]    poly_vtx         Polygon connectivity
 * \param [in]    mesh_vtx_coords  Coordinates of mesh vertices
 * \param [inout] mesh_vtx_coords  Mean value coordinates of points to locate
 *
 *
 */

void
PDM_mean_values_polygon_compute//--> REMOVE
(
 const int     n_pts,
 const double *pts,
 const int     n_vtx,
 const int    *poly_vtx,
 const double *mesh_vtx_coords,
 double       *mean_values
 )
{
  const double eps_base = 1e-10;

  printf("vtx coords =\n");
  double *vtx_coords = malloc (sizeof(double) * n_vtx * 3);
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    for (int idim = 0; idim < 3; idim++) {
      vtx_coords[3*ivtx + idim] = mesh_vtx_coords[3*poly_vtx[ivtx] + idim];
      printf("%f ", vtx_coords[3*ivtx + idim]);
    }
    printf("\n");
  }
  printf("\n\n");

  /* Projection onto average plane */
  double bary[3];
  PDM_polygon_compute_barycenter (n_vtx, vtx_coords, bary);
  printf("polygon barycenter = (%f, %f, %f)\n", bary[0], bary[1], bary[2]);

  double normal[3] = {0., 0., 1.};
  PDM_plane_normal (n_vtx, vtx_coords, normal);
  printf("polygon normal = (%f, %f, %f)\n", normal[0], normal[1], normal[2]);

  /* PDM_plane_projection2 (coo_point_dist, bary, n, coo_point_dist); */

  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    double *vco = vtx_coords + 3 * ivtx;
    PDM_plane_projection2 (vco, bary, normal, vco);
  }

  double poly_bounds[6] = {DBL_MAX, -DBL_MAX,
                           DBL_MAX, -DBL_MAX,
                           DBL_MAX, -DBL_MAX};
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    poly_bounds[0] = PDM_MIN (poly_bounds[0], vtx_coords[3*ivtx]);
    poly_bounds[1] = PDM_MAX (poly_bounds[1], vtx_coords[3*ivtx]);

    poly_bounds[2] = PDM_MIN (poly_bounds[2], vtx_coords[3*ivtx + 1]);
    poly_bounds[3] = PDM_MAX (poly_bounds[3], vtx_coords[3*ivtx + 1]);

    poly_bounds[4] = PDM_MIN (poly_bounds[4], vtx_coords[3*ivtx + 2]);
    poly_bounds[5] = PDM_MAX (poly_bounds[5], vtx_coords[3*ivtx + 2]);
  }
  printf("polygon bounds = [%f %f %f %f %f %f]\n\n\n",
         poly_bounds[0],
         poly_bounds[1],
         poly_bounds[2],
         poly_bounds[3],
         poly_bounds[4],
         poly_bounds[5]);
  // TO DO : EXPAND BOUNDS BY EPSILON


  /* Loop over points to locate */
  double _pt[3];
  double *s    = malloc (sizeof(double) * n_vtx * 3);
  double *dist = malloc (sizeof(double) * n_vtx);
  double *area = malloc (sizeof(double) * n_vtx);
  double *dotp = malloc (sizeof(double) * n_vtx);

  for (int ipt = 0; ipt < n_pts; ipt++) {

    for (int idim = 0; idim < 3; idim++) {
      _pt[idim] = pts[3*ipt + idim];
    }

    PDM_plane_projection2 (_pt, bary, normal, _pt);

    /* Make sure the point lies inside the polygon or on its boundary */
    if (PDM_polygon_point_in (_pt,
                              n_vtx,
                              vtx_coords,
                              poly_bounds,
                              normal) != 1) {

      printf("! point %d not in polygon (%f, %f, %f)\n",
             ipt,
             _pt[0],
             _pt[1],
             _pt[2]);

      double closest_point[3];
      double dist_min = DBL_MAX;

      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        double *v1 = vtx_coords + 3 * ivtx;
        double *v2 = vtx_coords + 3 * ((ivtx+1) % n_vtx);
        double closest[3];
        double t;

        double dist2 = PDM_line_distance (_pt,
                                          v1,
                                          v2,
                                          &t,
                                          closest);

        if (dist2 < dist_min) {
          dist_min = dist2;

          for (int idim = 0; idim < 3; idim++) {
            closest_point[idim] = closest[idim];
          }
        }
      }
      for (int idim = 0; idim < 3; idim++) {
        _pt[idim] = closest_point[idim];
      }


      printf("\t closest point = (%f, %f, %f), dist = %f\n",
             _pt[0],
             _pt[1],
             _pt[2],
             sqrt(dist_min));
    }


    /* Compute mean values */
    double min_dist = DBL_MAX;
    double *_mean_val = mean_values + n_vtx * ipt;
    int is_on_edge = 0;
    int is_vertex  = 0;

    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      int jvtx = (ivtx + 1) % n_vtx;

      double *vect = s + 3 * ivtx;
      for (int idim = 0; idim < 3; idim++) {
        vect[idim] = vtx_coords[3*jvtx + idim]- vtx_coords[3*ivtx + idim];
      }
      double l_edge = PDM_MODULE (vect);
      min_dist = PDM_MIN (l_edge, min_dist);
    }

    double eps = PDM_MAX (min_dist * eps_base, 1.e-30);

    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      double *vect = s + 3 * ivtx;
      for (int idim = 0; idim < 3; idim++) {
        vect[idim] = vtx_coords[3*ivtx + idim] - _pt[idim];
      }

      dist[ivtx] = PDM_MODULE (vect);
    }

    int current_vtx;
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      int jvtx = (ivtx + 1) % n_vtx;
      double *vect1 = s + 3 * ivtx;
      double *vect2 = s + 3 * jvtx;
      double crossp[3];

      dotp[ivtx] = PDM_DOT_PRODUCT (vect1, vect2);
      PDM_CROSS_PRODUCT(crossp, vect1, vect2);

      double sign = PDM_DOT_PRODUCT (crossp, normal);
      area[ivtx] = PDM_MODULE(crossp);

      if (sign < 0) {
        area[ivtx] = -area[ivtx];
      }

      if (dist[ivtx] <= eps) {

        is_vertex = 1;
        current_vtx = ivtx;
        break;

      } else if ((fabs(area[ivtx]) <= eps)  && (dotp[ivtx] < 0)) {

        is_on_edge = 1;
        current_vtx = ivtx;
        break;

      }
    }

    /* Current point is on a vertex of the polygon */
    if (is_vertex) {
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        _mean_val[ivtx] = 0.;
      }
      _mean_val[current_vtx] = 1.;
    }

    /* Current point is on an edge of the polygon */
    else if (is_on_edge) {
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        _mean_val[ivtx] = 0.;
      }

      int next_vtx = (current_vtx + 1) % n_vtx;
      double frac = 1. / (dist[current_vtx] + dist[next_vtx]);

      _mean_val[current_vtx] = dist[next_vtx]    * frac;
      _mean_val[next_vtx]    = dist[current_vtx] * frac;
    }

    /* General case */
    else {

      double sigma = 0.;
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        double coef = 0.;
        int prev_vtx = (ivtx - 1 + n_vtx) % n_vtx;
        int next_vtx = (ivtx + 1) % n_vtx;

        if (fabs(area[prev_vtx]) > eps) {
          coef += (dist[prev_vtx] - dotp[prev_vtx]/dist[ivtx]) / area[prev_vtx];
        }

        if (fabs(area[ivtx]) > eps) {
          coef += (dist[next_vtx] - dotp[ivtx]/dist[ivtx]) / area[ivtx];
        }

        sigma += coef;
        _mean_val[ivtx] = coef;
      }

      if (fabs(sigma) >= eps ) {

        double inv_sigma = 1. / sigma;
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          _mean_val[ivtx] *= inv_sigma;
        }

      } else {

        double abs_sigma = fabs(sigma);
        printf("Warning : Mise à NAN %f %f\n", abs_sigma,  eps);
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          _mean_val[ivtx] = NAN;
        }

      }

      /* Check result */
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        if (_mean_val[ivtx] != _mean_val[ivtx] ||
            _mean_val[ivtx] < 0. ||
            _mean_val[ivtx] > 1.) {

          double dist_min = DBL_MAX;
          int i_min = 0;
          double t_min;

          for (int i = 0; i < n_vtx; i++) {
            _mean_val[i] = 0.;
          }

          for (int i = 0; i < n_vtx; i++) {
            double *v1 = vtx_coords + 3 * i;
            double *v2 = vtx_coords + 3 * ((i+1) % n_vtx);
            double closest[3];
            double t;

            double dist2 = PDM_line_distance (_pt,
                                              v1,
                                              v2,
                                              &t,
                                              closest);

            if (dist2 < dist_min) {
              t_min = t;
              i_min = i;
            }
          }

          _mean_val[i_min] = 1 - t_min;
          _mean_val[(i_min + 1) % n_vtx] = t_min;

          break;
        }
      }

    }
  }

  free (vtx_coords);
  free (s);
  free (dist);
  free (area);
  free (dotp);
}






void
PDM_mean_values_polygon_compute2//--> REMOVE
(
 const int     n_pts,
 const double *pts,
 const int     n_vtx,
 const int    *poly_vtx,
 const double *mesh_vtx_coords,
 double       *mean_values
 )
{
  const double eps_base = 1e-10;

  double *vtx_xyz = malloc (sizeof(double) * n_vtx * 3);
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    for (int idim = 0; idim < 3; idim++) {
      vtx_xyz[3*ivtx + idim] = mesh_vtx_coords[3*poly_vtx[ivtx] + idim];
    }
  }

  double tangent_u[3], tangent_v[3], normal[3];
  const double *orig_xyz = mesh_vtx_coords + 3 * poly_vtx[0];

  PDM_polygon_orthonormal_basis (n_vtx,
                                 vtx_xyz,
                                 tangent_u,
                                 tangent_v,
                                 normal);

  double *vtx_uv = malloc (sizeof(double) * n_vtx * 2);
  PDM_polygon_compute_uv_coordinates (n_vtx,
                                      vtx_xyz,
                                      orig_xyz,
                                      tangent_u,
                                      tangent_v,
                                      vtx_uv);

  double uv_bounds[4] = {DBL_MAX, -DBL_MAX,
                         DBL_MAX, -DBL_MAX};

  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    uv_bounds[0] = PDM_MIN (uv_bounds[0], vtx_uv[2*ivtx]);
    uv_bounds[1] = PDM_MAX (uv_bounds[1], vtx_uv[2*ivtx]);

    uv_bounds[2] = PDM_MIN (uv_bounds[2], vtx_uv[2*ivtx + 1]);
    uv_bounds[3] = PDM_MAX (uv_bounds[3], vtx_uv[2*ivtx + 1]);
  }



  /* Loop over points to locate */
  double uv[2];
  double *s    = malloc (sizeof(double) * n_vtx * 2);
  double *dist = malloc (sizeof(double) * n_vtx);
  double *area = malloc (sizeof(double) * n_vtx);
  double *dotp = malloc (sizeof(double) * n_vtx);

  for (int ipt = 0; ipt < n_pts; ipt++) {
    PDM_polygon_compute_uv_coordinates (1,
                                        pts + 3*ipt,
                                        orig_xyz,
                                        tangent_u,
                                        tangent_v,
                                        uv);

    /* Make sure the point lies inside the polygon or on its boundary */
    if (PDM_polygon_point_in_2d (uv,
                                 n_vtx,
                                 vtx_uv,
                                 uv_bounds) != PDM_POLYGON_INSIDE) {
      double uv_closest[2];
      double dist_min = DBL_MAX;

      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        double t;
        double closest[2];
        double dist2 = PDM_line_distance_2d (uv,
                                             vtx_uv + 2*ivtx,
                                             vtx_uv + 2*((ivtx+1)%n_vtx),
                                             &t,
                                             closest);

        if (dist2 < dist_min) {
          dist_min = dist2;

          for (int idim = 0; idim < 2; idim++) {
            uv_closest[idim] = closest[idim];
          }
        }
      }

      for (int idim = 0; idim < 2; idim++) {
        uv[idim] = uv_closest[idim];
      }
    }


    /* Compute mean values */
    double min_dist = DBL_MAX;
    double *_mean_val = mean_values + n_vtx * ipt;
    int is_on_edge = 0;
    int is_vertex  = 0;

    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      int jvtx = (ivtx + 1) % n_vtx;

      double *vect = s + 2 * ivtx;
      for (int idim = 0; idim < 2; idim++) {
        vect[idim] = vtx_uv[2*jvtx + idim] - vtx_uv[2*ivtx + idim];
      }
      double l_edge = sqrt (PDM_DOT_PRODUCT_2D (vect, vect));
      min_dist = PDM_MIN (l_edge, min_dist);
    }

    double eps = PDM_MAX (min_dist * eps_base, 1.e-30);

    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      double *vect = s + 2 * ivtx;
      for (int idim = 0; idim < 2; idim++) {
        vect[idim] = vtx_uv[2*ivtx + idim] - uv[idim];
      }

      dist[ivtx] = sqrt (PDM_DOT_PRODUCT_2D (vect, vect));
    }

    int current_vtx;
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      int jvtx = (ivtx + 1) % n_vtx;
      double *vect1 = s + 2 * ivtx;
      double *vect2 = s + 2 * jvtx;

      dotp[ivtx] = PDM_DOT_PRODUCT_2D (vect1, vect2);
      area[ivtx] = vect1[0] * vect2[1] - vect1[1] * vect2[0];

      if (dist[ivtx] <= eps) {

        is_vertex = 1;
        current_vtx = ivtx;
        break;

      } else if ((fabs(area[ivtx]) <= eps)  && (dotp[ivtx] < 0)) {

        is_on_edge = 1;
        current_vtx = ivtx;
        break;

      }
    }

    /* Current point is on a vertex of the polygon */
    if (is_vertex) {
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        _mean_val[ivtx] = 0.;
      }
      _mean_val[current_vtx] = 1.;
    }

    /* Current point is on an edge of the polygon */
    else if (is_on_edge) {
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        _mean_val[ivtx] = 0.;
      }

      int next_vtx = (current_vtx + 1) % n_vtx;
      double frac = 1. / (dist[current_vtx] + dist[next_vtx]);

      _mean_val[current_vtx] = dist[next_vtx]    * frac;
      _mean_val[next_vtx]    = dist[current_vtx] * frac;
    }

    /* General case */
    else {

      double sum = 0.;
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        double coef = 0.;
        int prev_vtx = (ivtx - 1 + n_vtx) % n_vtx;
        int next_vtx = (ivtx + 1) % n_vtx;

        if (fabs(area[prev_vtx]) > eps) {
          coef += (dist[prev_vtx] - dotp[prev_vtx]/dist[ivtx]) / area[prev_vtx];
        }

        if (fabs(area[ivtx]) > eps) {
          coef += (dist[next_vtx] - dotp[ivtx]/dist[ivtx]) / area[ivtx];
        }

        sum += coef;
        _mean_val[ivtx] = coef;
      }

      if (fabs(sum) >= eps ) {

        sum = 1. / sum;
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          _mean_val[ivtx] *= sum;
        }

      } else {

        double abs_sigma = fabs(sum);
        printf("Warning : Mise à NAN %f %f\n", abs_sigma,  eps);
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          _mean_val[ivtx] = NAN;
        }

      }

      /* Check result */
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        if (_mean_val[ivtx] != _mean_val[ivtx] ||
            _mean_val[ivtx] < 0. ||
            _mean_val[ivtx] > 1.) {

          double dist_min = DBL_MAX;
          int i_min = 0;
          double t_min;

          for (int i = 0; i < n_vtx; i++) {
            _mean_val[i] = 0.;
          }

          for (int i = 0; i < n_vtx; i++) {
            double *v1 = vtx_uv + 2 * i;
            double *v2 = vtx_uv + 2 * ((i+1) % n_vtx);
            double closest[2];
            double t;

            double dist2 = PDM_line_distance_2d (uv,
                                                 v1,
                                                 v2,
                                                 &t,
                                                 closest);

            if (dist2 < dist_min) {
              t_min = t;
              i_min = i;
            }
          }

          _mean_val[i_min] = 1 - t_min;
          _mean_val[(i_min + 1) % n_vtx] = t_min;

          break;
        }
      }

    }
  }

  free (vtx_xyz);
  free (vtx_uv);
  free (s);
  free (dist);
  free (area);
  free (dotp);
}




void
PDM_mean_values_polygon_compute_2d//--> REMOVE
(
 const int     n_pts,
 const double *pts_xy,
 const int     n_vtx,
 const double *vtx_xy,
 double       *mean_value_coords
 )
{
  /*
    See "Mean value coordinates for arbitrary planar polygons", Kai Hormann, and Michael S. Floater. (2006)
   */
  const PDM_bool_t EXTEND_OUTSIDE = PDM_FALSE;//PDM_TRUE;
  const double eps_base = 1e-12;

  /* Compute polygon bounds */
  double bounds[4] = {DBL_MAX, -DBL_MAX,
                      DBL_MAX, -DBL_MAX};

  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    bounds[0] = PDM_MIN (bounds[0], vtx_xy[2*ivtx]);
    bounds[1] = PDM_MAX (bounds[1], vtx_xy[2*ivtx]);

    bounds[2] = PDM_MIN (bounds[2], vtx_xy[2*ivtx + 1]);
    bounds[3] = PDM_MAX (bounds[3], vtx_xy[2*ivtx + 1]);
  }

  double diag2 = (bounds[1] - bounds[0]) * (bounds[1] - bounds[0]) + (bounds[3] - bounds[2]) * (bounds[3] - bounds[2]);

  double eps_len2 = eps_base * eps_base * diag2;

  double *s = malloc (sizeof(double) * n_vtx * 2);
  double *r = malloc (sizeof(double) * n_vtx);
  double *A = malloc (sizeof(double) * n_vtx);
  double *D = malloc (sizeof(double) * n_vtx);

  for (int ipt = 0; ipt < n_pts; ipt++) {
    const double *xy = pts_xy + 2 * ipt;
    double *mvc = mean_value_coords + n_vtx * ipt;

    if (EXTEND_OUTSIDE == PDM_FALSE) {
      if (PDM_polygon_point_in_2d (xy,
                                   n_vtx,
                                   vtx_xy,
                                   bounds) != PDM_POLYGON_INSIDE) {
        /* Project point on polygon boundary */
        double dist_min = DBL_MAX;
        int imin;
        double tmin;

        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          double t;
          double closest[2];
          double dist2 = PDM_line_distance_2d (xy,
                                               vtx_xy + 2*ivtx,
                                               vtx_xy + 2*((ivtx+1)%n_vtx),
                                               &t,
                                               closest);

          if (dist2 < dist_min) {
            dist_min = dist2;
            imin = ivtx;
            tmin = t;
          }
        }

        mvc[imin]           = 1.0 - tmin;
        mvc[(imin+1)%n_vtx] = tmin;
        continue;
      }
    }

    PDM_bool_t special_case = PDM_FALSE;
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      double *vec = s + 2*ivtx;
      for (int idim = 0; idim < 2; idim++) {
        vec[idim] = vtx_xy[2*ivtx + idim] - xy[idim];
      }
      r[ivtx] = PDM_DOT_PRODUCT_2D (vec, vec);

      if (r[ivtx] < eps_len2) {
        /* Point coincident with vertex */
        for (int i = 0; i < n_vtx; i++) {
          mvc[i] = 0.0;
        }
        mvc[ivtx] = 1.0;

        special_case = PDM_TRUE;
        break;
      }

      r[ivtx] = sqrt (r[ivtx]);
    }

    if (special_case == PDM_TRUE) {
      continue;
    }

    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      int jvtx = (ivtx + 1) % n_vtx;
      double *veci = s + 2*ivtx;
      double *vecj = s + 2*jvtx;

      A[ivtx] = veci[0] * vecj[1] - veci[1] * vecj[0];
      D[ivtx] = PDM_DOT_PRODUCT_2D (veci, vecj);

      if (fabs(A[ivtx]) < eps_base && D[ivtx] < 0) {
        /* Point on edge */
        double vec[2];
        vec[0] = vtx_xy[2*jvtx] - vtx_xy[2*ivtx];
        vec[1] = vtx_xy[2*jvtx+1] - vtx_xy[2*ivtx+1];

        double lvec2 = PDM_DOT_PRODUCT_2D (vec, vec);
        double t;
        if (lvec2 < eps_len2) {
          /* Degenerate edge */
          t = 0.0;
        } else {
          t = - PDM_DOT_PRODUCT_2D (veci, vec) / lvec2;
          //t = r[ivtx] / sqrt(lvec2);

          if (t < -eps_base || t > 1.0 + eps_base) {
            /* Error */
            /* ... */
          }
        }

        //printf("Point on edge: t = %f, 1 - t = %f\n", t, 1 - t);

        for (int i = 0; i < n_vtx; i++) {
          mvc[i] = 0.0;
        }
        mvc[ivtx] = 1.0 - t;
        mvc[jvtx] = t;

        special_case = PDM_TRUE;
        break;

      }
    }

    if (special_case == PDM_TRUE) {
      continue;
    }

    /* General case (point strictly inside polygon) */
    double sum_w = 0.0;
    for (int i = 0; i < n_vtx; i++) {
      int ip = (i + 1) % n_vtx;
      int im = (i - 1 + n_vtx) % n_vtx;

      double w = 0.0;
      if (fabs(A[i]) > eps_base) {
        w += (r[ip] - D[i] / r[i]) / A[i];
      }

      if (fabs(A[im]) > eps_base) {
        w += (r[im] - D[im] / r[i]) / A[im];
      }

      sum_w += w;
      mvc[i] = w;
    }

    if (fabs(sum_w) > eps_base) {
      sum_w = 1.0 / sum_w;
      for (int i = 0; i < n_vtx; i++) {
        mvc[i] *= sum_w;
      }
    }

  }
}




void
PDM_mean_values_polygon_compute3//--> REMOVE
(
 const int     n_pts,
 const double *pts_xyz,
 const int     n_vtx,
 const double *vtx_xyz,
 double       *mean_value_coords
 )
{
  double tangent_u[3], tangent_v[3], normal[3];
  const double *orig_xyz = vtx_xyz;

  PDM_polygon_orthonormal_basis (n_vtx,
                                 vtx_xyz,
                                 tangent_u,
                                 tangent_v,
                                 normal);

  double *vtx_uv = malloc (sizeof(double) * n_vtx * 2);
  PDM_polygon_compute_uv_coordinates (n_vtx,
                                      vtx_xyz,
                                      orig_xyz,
                                      tangent_u,
                                      tangent_v,
                                      vtx_uv);

  double *pts_uv = malloc (sizeof(double) * n_pts * 2);
  PDM_polygon_compute_uv_coordinates (n_pts,
                                      pts_xyz,
                                      orig_xyz,
                                      tangent_u,
                                      tangent_v,
                                      pts_uv);

#if 0
  PDM_mean_values_polygon_compute_2d (n_pts,
                                      pts_uv,
                                      n_vtx,
                                      vtx_uv,
                                      mean_value_coords);
#else
  PDM_mean_value_coordinates_polygon_2d (n_vtx,
                                         vtx_uv,
                                         n_pts,
                                         pts_uv,
                                         mean_value_coords);
#endif
}









void
PDM_mean_values_polyhedron
(
 const int     n_pts,
 const double  pts_xyz[],
 const int     n_vtx,
 const double  vtx_xyz[],
 const int     n_faces,
 const int     face_vtx_idx[],
 const int     face_vtx[],
 double       *mean_value_coords
 )
{
  double *vtx_xyz_o = malloc (sizeof(double) * n_vtx * 3);
  double *mag_vtx_xyz_o = malloc (sizeof(double) * n_vtx);
  double *vtx_s = malloc (sizeof(double) * n_vtx * 3);
  double *nor = malloc (sizeof(double) * n_vtx * 3);

  double uv_origin[2] = {0.0, 0.0};

  const double epsilon = 1e-10; // Adapt to polyhedron scale?
  const double epsilon2 = epsilon * epsilon;

  for (int ipt = 0; ipt < n_pts; ipt++) {
    const double *p_xyz = pts_xyz + 3 * ipt;
    double *mvc = mean_value_coords + n_vtx * ipt;
    for (int i = 0; i < n_vtx; i++) {
      mvc[i] = 0.0;
    }

    /*printf("\n");
    for (int k = 0; k < n_vtx; k++) {
      printf(" %f", mvc[k]);
    }
    printf("\n");*/

    PDM_bool_t special_case = PDM_FALSE;

    /* Offset polyhedron's vertices */
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      double *vec = vtx_xyz_o + 3 * ivtx;
      mag_vtx_xyz_o[ivtx] = 0.;

      for (int idim = 0; idim < 3; idim++) {
        vec[idim] = vtx_xyz[3*ivtx + idim] - p_xyz[idim];
        mag_vtx_xyz_o[ivtx] += vec[idim] * vec[idim];
      }

      if (mag_vtx_xyz_o[ivtx] < epsilon2) {
        /* Point coincident with vertex ivtx */
        mvc[ivtx] = 1.0;
        special_case = PDM_TRUE;
        printf("Point #%d on vertex %d\n", ipt, ivtx);

        printf("Point #%d special case, mvc = (", ipt);
        for (int k = 0; k < n_vtx; k++) {
          printf(" %f", mvc[k]);
        }
        printf(")\n\n");
        break; /* exit loop over vertices */

      } else {
        mag_vtx_xyz_o[ivtx] = 1.0 / sqrt(mag_vtx_xyz_o[ivtx]);
      }

    }

    if (special_case == PDM_TRUE) continue; /* continue loop over points */


    /* Project offset vertices onto unit sphere */
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      for (int idim = 0; idim < 3; idim++) {
        vtx_s[3*ivtx + idim] = vtx_xyz_o[3*ivtx + idim] * mag_vtx_xyz_o[ivtx];
      }
    }


    /* Loop over faces */
    for (int iface = 0; iface < n_faces; iface++) {
      const int *f = face_vtx + face_vtx_idx[iface];
      int n_vtx_f = face_vtx_idx[iface+1] - face_vtx_idx[iface];

      /* Compute face vector vf */
      double vf[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < n_vtx_f; i++) {
        int ip = (i + 1) % n_vtx_f;

        double *n_i = nor + 3*i;
        double *v_i = vtx_s + 3*f[i];
        double *v_ip = vtx_s + 3*f[ip];
        PDM_CROSS_PRODUCT (n_i, v_i, v_ip);

        double mag_nor = PDM_DOT_PRODUCT (n_i, n_i);

        if (mag_nor < epsilon2) {
          /* Point on edge (f[i], f[ip]) */
          double vec[3];
          double lvec2 = 0.0;
          for (int idim = 0; idim < 3; idim++) {
            vec[idim] = vtx_xyz[3*f[ip] + idim] - vtx_xyz[3*f[i] + idim];
            lvec2 += vec[idim] * vec[idim];
          }

          double t = 0.0;
          if (lvec2 > epsilon2) {
            double *v = vtx_xyz_o + 3*f[i];
            t = -PDM_DOT_PRODUCT (v, vec) / lvec2;
          }
          for (int k = 0; k < n_vtx; k++) {
            mvc[k] = 0.0;
          }

          mvc[f[i]] = 1.0 - t;
          mvc[f[ip]] = t;
          special_case = PDM_TRUE;
          break; /* exit loop over face vertices */

        } else {

          for (int idim = 0; idim < 3; idim++) {
            vf[idim] += 0.5 * n_i[idim];
          }
        }

      } /* Loop over face vertices */

      if (special_case == PDM_TRUE) break; /* exit loop over faces */


      /* Compute generalized barycentric coords of vf in current face */
      /*    Projection onto tangent plane */
      double mag_vf = PDM_DOT_PRODUCT (vf, vf);
      double p_s[3];

      if (mag_vf < epsilon2) {
        /* Error? */
        printf("!!! |vf| << 1\n");
        p_s[0] = 0.0;
        p_s[1] = 0.0;
        p_s[2] = 0.0;

      } else {
        mag_vf = 1.0 / sqrt(mag_vf);
        for (int idim = 0; idim < 3; idim++) {
          p_s[idim] = vf[idim] * mag_vf;
        }
      }

      double *vtx_f_p  = malloc (sizeof(double) * n_vtx_f * 3);
      double *mvc_f    = malloc (sizeof(double) * n_vtx_f);

      for (int i = 0; i < n_vtx_f; i++) {
        int j = f[i];

        double *v = vtx_s + 3*j;
        double denom = PDM_DOT_PRODUCT (v, p_s);

        if (fabs(denom) < epsilon) {
          /* Point on face */
          for (int k = 0; k < n_vtx_f; k++) {
            int ivtx = f[k];
            for (int idim = 0; idim < 3; idim++) {
              vtx_f_p[3*k + idim] = vtx_xyz[3*ivtx + idim];
            }
          }
          PDM_mean_values_polygon_compute3 (1,
                                            p_xyz,
                                            n_vtx_f,
                                            vtx_f_p,
                                            mvc_f);

          for (int k = 0; k < n_vtx; k++) {
            mvc[k] = 0.0;
          }
          for (int k = 0; k < n_vtx_f; k++) {
            int ivtx = f[k];
            mvc[ivtx] = mvc_f[k];
          }

          special_case = PDM_TRUE;
          break;
        }
        denom = 1.0 / denom;

        for (int idim = 0; idim < 3; idim++) {
          vtx_f_p[3*i + idim] = v[idim] * denom;
        }
      }

      if (special_case == PDM_TRUE) {
        free (vtx_f_p);
        free (mvc_f);
        break;
      }

      /*    Compute uv coords in tangent plane */
      double *v0 = vtx_f_p;
      double *v1 = vtx_f_p + 3;

      double tangent_u[3], tangent_v[3];
      double mag2 = 0;
      for (int idim = 0; idim < 3; idim++) {
        tangent_u[idim] = v1[idim] - v0[idim];
        mag2 += tangent_u[idim] * tangent_u[idim];
      }
      mag2 = 1.0 / sqrt(mag2);

      for (int idim = 0; idim < 3; idim++) {
        tangent_u[idim] *= mag2;
      }

      PDM_CROSS_PRODUCT (tangent_v, p_s, tangent_u);

      double *vtx_f_uv = malloc (sizeof(double) * n_vtx_f * 2);
      PDM_polygon_compute_uv_coordinates (n_vtx_f,
                                          vtx_f_p,
                                          p_s,
                                          tangent_u,
                                          tangent_v,
                                          vtx_f_uv);

      /*    Compute planar mean value coords */
      PDM_mean_values_polygon_compute_2d (1,
                                          uv_origin,
                                          n_vtx_f,
                                          vtx_f_uv,
                                          mvc_f);
      free (vtx_f_uv);

      /*    Compute generalized mean value coords and transfer weights to vertices */
      for (int i = 0; i < n_vtx_f; i++) {
        int j = f[i];
        double *vsj = vtx_s + 3*j;
        double *vpi = vtx_f_p + 3*i;

        mvc_f[i] *= PDM_DOT_PRODUCT (vsj, vpi) * mag_vf * mag_vtx_xyz_o[j];
        mvc[j] += mvc_f[i];
      }
      free (mvc_f);
      free (vtx_f_p);

    } /* Loop over faces */


    if (special_case == PDM_FALSE) {
      /* Mean value coords in polyhedron */
      double sum_w = 0.0;
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        sum_w += mvc[ivtx];
      }

      if (fabs(sum_w) < epsilon) {
        /* Error? */
        printf("Point #%d (%f %f %f) : |Sum(w)| << 1\n", ipt, p_xyz[0], p_xyz[1], p_xyz[2]);
      } else {
        sum_w = 1.0 / sum_w;

        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          mvc[ivtx] *= sum_w;
        }
      }

    }

  } /* Loop over points */

  free (vtx_xyz_o);
  free (mag_vtx_xyz_o);
  free (vtx_s);
  free (nor);
}









void
PDM_mean_value_coordinates_polygon_2d
(
 const int    n_vtx,
 const double vtx_coord[],
 const int    n_pts,
 const double pts_coord[],
 double       mean_value_coord[]
)
{
  /*
    See "Mean value coordinates for arbitrary planar polygons", Kai Hormann, and Michael S. Floater. (2006)
  */
  const double eps_base = 1e-12;

  int ipt, ivtx, jvtx, idim;

  /* Compute polygon bounds */
  double bounds[4] = {DBL_MAX, -DBL_MAX,
                      DBL_MAX, -DBL_MAX};
  double char_length = 0.;
  for (ivtx = 0; ivtx < n_vtx; ivtx++) {
    for (idim = 0; idim < 2; idim++) {
      bounds[2*idim]   = PDM_MIN (bounds[2*idim],   vtx_coord[2*ivtx+idim]);
      bounds[2*idim+1] = PDM_MAX (bounds[2*idim+1], vtx_coord[2*ivtx+idim]);

      char_length = PDM_MAX (char_length, bounds[2*idim+1] - bounds[2*idim+1]);
    }
  }

  double eps2 = eps_base * char_length;
  eps2 *= eps2;

  double *s = malloc (sizeof(double) * n_vtx * 2);
  double *r = malloc (sizeof(double) * n_vtx);
  double *A = malloc (sizeof(double) * n_vtx);
  double *D = malloc (sizeof(double) * n_vtx);

  /* Loop on points */
  for (ipt = 0; ipt < n_pts; ipt++) {
    const double *_pt = pts_coord + 2 * ipt;
    double       *_bc = mean_value_coord + n_vtx * ipt;

    /* If current point is outside the polygon,
       consider its projection on the polygon boundary */
    if (PDM_polygon_point_in2d (_pt,
                                n_vtx,
                                vtx_coord,
                                char_length,
                                bounds) != PDM_POLYGON_INSIDE) {
      double dist2_min = DBL_MAX;
      int i_min;
      double t_min;

      for (ivtx = 0; ivtx < n_vtx; ivtx++) {
        jvtx = (ivtx + 1) % n_vtx;
        double t;
        double closest[2];
        double dist2 = PDM_line_distance_2d (_pt,
                                             vtx_coord + 2*ivtx,
                                             vtx_coord + 2*jvtx,
                                             &t,
                                             closest);

        if (dist2 < dist2_min) {
          dist2_min = dist2;
          i_min = ivtx;

          if (t < 0.) {
            t_min = 0.;
          } else if (t > 1.) {
            t_min = 1.;
          } else  {
            t_min = t;
          }
        }
      }

      for (int i = 0; i < n_vtx; i++) {
        _bc[i] = 0.;
      }
      _bc[i_min]           = 1.0 - t_min;
      _bc[(i_min+1)%n_vtx] = t_min;
      continue;
    }

    PDM_bool_t special_case = PDM_FALSE;
    for (ivtx = 0; ivtx < n_vtx; ivtx++) {
      double *vec = s + 2*ivtx;
      for (idim = 0; idim < 2; idim++) {
        vec[idim] = vtx_coord[2*ivtx + idim] - _pt[idim];
      }
      r[ivtx] = PDM_DOT_PRODUCT_2D (vec, vec);

      if (r[ivtx] < eps2) {
        /* Point coincident with vertex */
        for (int i = 0; i < n_vtx; i++) {
          _bc[i] = 0.;
        }
        _bc[ivtx] = 1.;

        special_case = PDM_TRUE;
        break;
      }

      else {
        r[ivtx] = sqrt (r[ivtx]);
      }
    }


    if (special_case == PDM_TRUE) {
      continue;
    }

    for (ivtx = 0; ivtx < n_vtx; ivtx++) {
      jvtx = (ivtx + 1) % n_vtx;
      double *veci = s + 2*ivtx;
      double *vecj = s + 2*jvtx;

      A[ivtx] = veci[0] * vecj[1] - veci[1] * vecj[0];
      D[ivtx] = PDM_DOT_PRODUCT_2D (veci, vecj);

      /* Point on line supporting edge */
      if (fabs(A[ivtx]) < eps_base && D[ivtx] < 0) {
        for (int i = 0; i < n_vtx; i++) {
          _bc[i] = 0.;
        }

#if 1
        double inv_denom = 1. / (r[ivtx] + r[jvtx]);
        _bc[ivtx] = r[jvtx] * inv_denom;
        _bc[jvtx] = r[ivtx] * inv_denom;
#else
        double vec[2];
        vec[0] = vtx_coord[2*jvtx] - vtx_coord[2*ivtx];
        vec[1] = vtx_coord[2*jvtx+1] - vtx_coord[2*ivtx+1];

        double lvec2 = PDM_DOT_PRODUCT_2D (vec, vec);
        double t;
        if (lvec2 < eps2) {
          /* Degenerate edge */
          t = 0.0;
        } else {
          t = - PDM_DOT_PRODUCT_2D (veci, vec) / lvec2;
        }

        _bc[ivtx] = 1. - t;
        _bc[jvtx] =      t;
#endif

        special_case = PDM_TRUE;
        break;
      }
    }

    if (special_case == PDM_TRUE) {
      continue;
    }

    /* General case (point strictly inside polygon) */
    double sum_w = 0.0;
    for (int i = 0; i < n_vtx; i++) {
      int ip = (i + 1) % n_vtx;
      int im = (i - 1 + n_vtx) % n_vtx;

      double w = 0.0;
      if (fabs(A[i]) > eps_base) {
        w += (r[ip] - D[i] / r[i]) / A[i];
      }

      if (fabs(A[im]) > eps_base) {
        w += (r[im] - D[im] / r[i]) / A[im];
      }

      sum_w += w;
      _bc[i] = w;
    }

    if (fabs(sum_w) > eps_base) {
      sum_w = 1.0 / sum_w;
      for (int i = 0; i < n_vtx; i++) {
        _bc[i] *= sum_w;
      }
    }

    else {
      printf("!!! sum_w = %g\n", sum_w);
    }
  }
}




void
PDM_mean_value_coordinates_polygon_3d
(
 const int    n_vtx,
 const double vtx_coord[],
 const int    n_pts,
 const double pts_coord[],
 double       mean_value_coord[]
)
{

  double *vtx_uv = malloc (sizeof(double) * n_vtx * 2);
  double *pts_uv = malloc (sizeof(double) * n_pts * 2);

  PDM_polygon_3d_to_2d (n_vtx,
                        vtx_coord,
                        vtx_uv,
                        n_pts,
                        pts_coord,
                        pts_uv,
                        NULL);

  PDM_mean_value_coordinates_polygon_2d (n_vtx,
                                         vtx_uv,
                                         n_pts,
                                         pts_uv,
                                         mean_value_coord);

  free (vtx_uv);
  free (pts_uv);
}




void
PDM_mean_value_coordinates_polyhedron
(
 const int         n_vtx,
 const double      vtx_coord[],
 const PDM_l_num_t n_face,
 const PDM_l_num_t face_vtx_idx[],
 const PDM_l_num_t face_vtx[],
 const int         face_orientation[],
 const double      pt_coord[],
 double            mean_value_coord[]
 )
{
  const int DEBUG = 0;

  if (DEBUG) {
    printf("\n\n\n--- PDM_mean_value_coordinates_polyhedron ---\n");
    printf("p_coord = %f %f %f\n", pt_coord[0], pt_coord[1], pt_coord[2]);
  }

  const double eps = 1.e-6;
  const double eps2 = eps * eps;

  double *v = malloc (sizeof(double) * n_vtx * 3);
  double *inv_lv = malloc (sizeof(double) * n_vtx);
  double *e = malloc (sizeof(double) * n_vtx * 3);

  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    mean_value_coord[ivtx] = 0.;
  }

  /*
   *  Offset vertices and project on unit sphere centered at point to locate
   */
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    double *_v = v + 3*ivtx;
    double lv2 = 0.;
    for (int idim = 0; idim < 3; idim++) {
      _v[idim] = vtx_coord[3*ivtx + idim] - pt_coord[idim];
      lv2 += _v[idim] * _v[idim];
    }

    /* Point coincident with vertex */
    if (lv2 < eps2) {
      mean_value_coord[ivtx] = 1.;

      free (v);
      free (inv_lv);
      free (e);
      return;
    }

    inv_lv[ivtx] = 1. / sqrt(lv2);

    double *_e = e + 3*ivtx;
    for (int idim = 0; idim < 3; idim++) {
      _e[idim] = _v[idim] * inv_lv[ivtx];
    }
  } // End of loop on vertices



  /*
   *  Prepare face triangulation
   */
  /* Count max nb of vertices per face */
  PDM_l_num_t n_vtx_face, n_vtx_face_max = 0;
  for (int iface = 0; iface < n_face; iface++) {
    n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];
    if (n_vtx_face > n_vtx_face_max) {
      n_vtx_face_max = n_vtx_face;
    }
  }

  /*
   *  Loop on faces
   */
  double m[3][3], b[3], denom[3];
  PDM_l_num_t n_tri;
  PDM_l_num_t _tri_vtx[3];
  PDM_l_num_t *tri_vtx = malloc (sizeof(PDM_l_num_t) * (n_vtx_face_max - 2)*3);
  PDM_triangulate_state_t *state = PDM_triangulate_state_create (n_vtx_face_max);

  for (int iface = 0; iface < n_face; iface++) {

    if (DEBUG) {
      printf("iface = %d\n", iface);
    }
    const PDM_l_num_t *_face_vtx = face_vtx + face_vtx_idx[iface];

    /*
     *  Triangulate current face
     */
    n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

    /* Triangular face */
    if (n_vtx_face == 3) {
      n_tri = 1;
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        tri_vtx[ivtx] = _face_vtx[ivtx];
      }
    }

    /* Quadrilateral face */
    else if (n_vtx_face == 4) {
      n_tri = PDM_triangulate_quadrangle (3,
                                          vtx_coord,
                                          NULL,
                                          _face_vtx,
                                          tri_vtx);
    }

    /* Polygonal face */
    else {
      n_tri = PDM_triangulate_polygon(3,
                                      n_vtx_face,
                                      vtx_coord,
                                      NULL,
                                      _face_vtx,
                                      PDM_TRIANGULATE_MESH_DEF,
                                      tri_vtx,
                                      state);
    }


    /* Loop on triangles */
    for (int itri = 0; itri < n_tri; itri++) {
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        _tri_vtx[ivtx] = tri_vtx[3*itri + ivtx] - 1;
      }

      if (face_orientation[iface] < 0) {
        PDM_l_num_t tmp = _tri_vtx[0];
        _tri_vtx[0] = _tri_vtx[2];
        _tri_vtx[2] = tmp;
      }

      if (DEBUG) {
        printf("   itri = %d, _tri_vtx = [%d %d %d]\n",
               itri, _tri_vtx[0], _tri_vtx[1], _tri_vtx[2]);
      }

      for (int j = 0; j < 3; j++) {
        int i = _tri_vtx[j];
        int ip = _tri_vtx[(j+1)%3];

        PDM_CROSS_PRODUCT (m[j], (e + 3*i), (e + 3*ip));
        double lm = PDM_DOT_PRODUCT (m[j], m[j]);

        if (lm < eps2) {
          printf("!!! (%f %f %f) face %d, lm = %g\n",
                 pt_coord[0], pt_coord[1], pt_coord[2],
                 iface,
                 sqrt(lm));
        }

        lm = sqrt(lm);
        b[j] = asin(lm);

        lm = 1. / lm;
        for (int idim = 0; idim < 3; idim++) {
          m[j][idim] *= lm;
        }
      }

      for (int i = 0; i < 3; i++) {
        denom[i] = PDM_DOT_PRODUCT ((e + 3*_tri_vtx[i]), m[(i+1)%3]);

        if (fabs(denom[i]) < eps) {
          /* Point located on current face(?) */
          n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];
          double *face_coord = malloc (sizeof(double) * n_vtx_face * 3);
          double *mean_value_coord_face = malloc (sizeof(double) * n_vtx_face);

          for (int j = 0; j < n_vtx_face; j++) {
            int id_vtx = face_vtx[face_vtx_idx[iface] + j] - 1;
            for (int idim = 0; idim < 3; idim++) {
              face_coord[3*j + idim] = vtx_coord[3*id_vtx + idim];
            }
          }

          PDM_mean_value_coordinates_polygon_3d (n_vtx_face,
                                                 face_coord,
                                                 1,
                                                 pt_coord,
                                                 mean_value_coord_face);

          for (int j = 0; j < n_vtx; j++) {
            mean_value_coord[j] = 0.;
          }

          for (int j = 0; j < n_vtx_face; j++) {
            int id_vtx = face_vtx[face_vtx_idx[iface] + j] - 1;
            mean_value_coord[id_vtx] = mean_value_coord_face[j];
          }

          free (face_coord);
          free (mean_value_coord_face);
          free (v);
          free (inv_lv);
          free (e);
          free (tri_vtx);
          return;
        }

        denom[i] = 0.5 / denom[i];

      }

      for (int i = 0; i < 3; i++) {
        int ivtx = _tri_vtx[i];
        int ip = (i+1)%3;

        for (int j = 0; j < 3; j++) {
          mean_value_coord[ivtx] += b[j] * denom[i] * PDM_DOT_PRODUCT (m[j], m[ip]);
        }
      }
    } // End of loop on triangles

  } // End of loop on faces


  /* Normalize */
  double sum = 0.;
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    mean_value_coord[ivtx] *= inv_lv[ivtx];
    sum += mean_value_coord[ivtx];
  }

  if (fabs(sum) > 1.e-12) {
    sum = 1. / sum;
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      mean_value_coord[ivtx] *= sum;
    }
  }

  free (v);
  free (inv_lv);
  free (e);
  free (tri_vtx);

}





void
PDM_mean_value_coordinates_polyhedron2
(
 const int         n_vtx,
 const double      vtx_coord[],
 const PDM_l_num_t n_face,
 const PDM_l_num_t face_vtx_idx[],
 const PDM_l_num_t face_vtx[],
 const int         face_orientation[],
 const double      pt_coord[],
 double            mean_value_coord[]
 )
{
  const double eps = 1.e-6;
  const double eps2 = eps * eps;

  const double eps_face = 1e-9;

  double *v = malloc (sizeof(double) * n_vtx * 3);
  double *inv_lv = malloc (sizeof(double) * n_vtx);
  double *e = malloc (sizeof(double) * n_vtx * 3);

  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    mean_value_coord[ivtx] = 0.;
  }

  /*
   *  Offset vertices and project on unit sphere centered at point to locate
   */
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    double *_v = v + 3*ivtx;
    double lv2 = 0.;
    for (int idim = 0; idim < 3; idim++) {
      _v[idim] = vtx_coord[3*ivtx + idim] - pt_coord[idim];
      lv2 += _v[idim] * _v[idim];
    }

    /* Point coincident with vertex */
    if (lv2 < eps2) {
      mean_value_coord[ivtx] = 1.;

      free (v);
      free (inv_lv);
      free (e);
      return;
    }

    inv_lv[ivtx] = 1. / sqrt(lv2);

    double *_e = e + 3*ivtx;
    for (int idim = 0; idim < 3; idim++) {
      _e[idim] = _v[idim] * inv_lv[ivtx];
    }
  } // End of loop on vertices



  /*
   *  Prepare face triangulation
   */
  /* Count max nb of vertices per face */
  PDM_l_num_t n_vtx_face, n_vtx_face_max = 0;
  for (int iface = 0; iface < n_face; iface++) {
    n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];
    if (n_vtx_face > n_vtx_face_max) {
      n_vtx_face_max = n_vtx_face;
    }
  }

  /*
   *  Loop on faces
   */
  double m[3][3], b[3];
  PDM_l_num_t n_tri;
  PDM_l_num_t _tri_vtx[3];
  PDM_l_num_t *tri_vtx = malloc (sizeof(PDM_l_num_t) * (n_vtx_face_max - 2)*3);
  PDM_triangulate_state_t *state = PDM_triangulate_state_create (n_vtx_face_max);

  for (int iface = 0; iface < n_face; iface++) {

    const PDM_l_num_t *_face_vtx = face_vtx + face_vtx_idx[iface];

    /*
     *  Triangulate current face
     */
    n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

    /* Triangular face */
    if (n_vtx_face == 3) {
      n_tri = 1;
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        tri_vtx[ivtx] = _face_vtx[ivtx];
      }
    }

    /* Quadrilateral face */
    else if (n_vtx_face == 4) {
      n_tri = PDM_triangulate_quadrangle (3,
                                          vtx_coord,
                                          NULL,
                                          _face_vtx,
                                          tri_vtx);
    }

    /* Polygonal face */
    else {
      n_tri = PDM_triangulate_polygon(3,
                                      n_vtx_face,
                                      vtx_coord,
                                      NULL,
                                      _face_vtx,
                                      PDM_TRIANGULATE_MESH_DEF,
                                      tri_vtx,
                                      state);
    }


    /* Loop on triangles */
    for (int itri = 0; itri < n_tri; itri++) {
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        _tri_vtx[ivtx] = tri_vtx[ivtx] - 1;
      }

      if (face_orientation[iface] < 0) {
        PDM_l_num_t tmp = _tri_vtx[0];
        _tri_vtx[0] = _tri_vtx[2];
        _tri_vtx[2] = tmp;
      }

      // Check triangle area...

      for (int isom = 0; isom < 3; isom++) {
        int isuiv = _tri_vtx[(isom+1)%3];
        int iprec = _tri_vtx[(isom+2)%3];

        double prod_scal = e[3*iprec    ] * e[3*isuiv    ]
          + e[3*iprec + 1] * e[3*isuiv + 1]
          + e[3*iprec + 2] * e[3*isuiv + 2];

        b[isom] = acos(prod_scal);

        m[isom][0] =  e[3*iprec + 1] * e[3*isuiv + 2]
          - e[3*iprec + 2] * e[3*isuiv + 1];
        m[isom][1] =  e[3*iprec + 2] * e[3*isuiv    ]
          - e[3*iprec    ] * e[3*isuiv + 2];
        m[isom][2] =  e[3*iprec    ] * e[3*isuiv + 1]
          - e[3*iprec + 1] * e[3*isuiv    ];

        double mod = sqrt(m[isom][0] * m[isom][0]
                          + m[isom][1] * m[isom][1]
                          + m[isom][2] * m[isom][2]);

        if (mod <  eps_face) {
          m[isom][0] = 0.;
          m[isom][1] = 0.;
          m[isom][2] = 0.;
        }

        else {

          m[isom][0] /= mod;
          m[isom][1] /= mod;
          m[isom][2] /= mod;
        }
      }

      for (int isom = 0; isom < 3; isom++) {

        double ps_nij_njk; //a ameliorer
        double ps_nki_njk; //a ameliorer
        double ps_ei_njk;  //a ameliorer

        const int iprec = (isom + 2) % 3;
        const int isuiv = (isom + 1) % 3;

        ps_nij_njk = m[isom][0] * m[isuiv][0]
        + m[isom][1] * m[isuiv][1]
        + m[isom][2] * m[isuiv][2];

        ps_nki_njk = m[isom][0] * m[iprec][0]
          + m[isom][1] * m[iprec][1]
          + m[isom][2] * m[iprec][2];

        // ps_ei_njk --> sur la face


        const int ivertex_tri = _tri_vtx[isom];
        ps_ei_njk = e[3*ivertex_tri    ] * m[isom][0]
          + e[3*ivertex_tri + 1] * m[isom][1]
          + e[3*ivertex_tri + 2] * m[isom][2];

        // vérifier ps_ei_njk

        if (fabs(ps_ei_njk) >  eps_face) {
          mean_value_coord[ivertex_tri] +=
            (b[isom] + b[isuiv] * ps_nij_njk + b[iprec] * ps_nki_njk)
            / (2 * ps_ei_njk);
        }

      } // Loop on vertices

    } // End of loop on triangles

  } // End of loop on faces


  /* Normalize */
  double sum = 0.;
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    mean_value_coord[ivtx] *= inv_lv[ivtx];
    sum += mean_value_coord[ivtx];
  }

  if (fabs(sum) > 1.e-12) {
    sum = 1. / sum;
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      mean_value_coord[ivtx] *= sum;
    }
  }

  free (v);
  free (inv_lv);
  free (e);
  free (tri_vtx);

}


static inline double
_determinant_3x3
(
 const double a[3],
 const double b[3],
 const double c[3]
 )
{
  return a[0] * (b[1]*c[2] - b[2]*c[1])
    +    a[1] * (b[2]*c[0] - b[0]*c[2])
    +    a[2] * (b[0]*c[1] - b[1]*c[0]);
}


/**
 * See "Mean value coordinates for closed triangular meshes", T. Ju et al. (2005)
 *
 *
 **/
void
PDM_mean_value_coordinates_polyhedron3
(
 const int         n_vtx,
 const double      vtx_coord[],
 const PDM_l_num_t n_face,
 const PDM_l_num_t face_vtx_idx[],
 const PDM_l_num_t face_vtx[],
 const int         face_orientation[],
 const double      pt_coord[],
 double            mean_value_coord[]
 )
{
  const int DEBUG = 0;//(pt_coord[0] < 0);
  if (DEBUG) {
    printf("\n\n-- PDM_mean_value_coordinates_polyhedron3 --\n");
    printf("pt_coord = %f %f %f\n", pt_coord[0], pt_coord[1], pt_coord[2]);
  }
  const int LOCATE_ON_TRIANGLES = 0;

  const double eps = 1.e-9;
  const double eps2 = eps * eps;

  double *u = malloc (sizeof(double) * n_vtx * 3);
  double *d = malloc (sizeof(double) * n_vtx);

  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    mean_value_coord[ivtx] = 0.;
  }

  /*
   *  Offset vertices and project on unit sphere centered at point to locate
   */
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    double *_u = u + 3*ivtx;
    for (int idim = 0; idim < 3; idim++) {
      _u[idim] = vtx_coord[3*ivtx + idim] - pt_coord[idim];
    }
    double uu = PDM_DOT_PRODUCT (_u, _u);

    /* Point coincident with vertex */
    if (uu < eps2) {
      mean_value_coord[ivtx] = 1.;

      free (u);
      free (d);
      return;
    }

    d[ivtx] = sqrt(uu);

    for (int idim = 0; idim < 3; idim++) {
      _u[idim] = _u[idim] / d[ivtx];
    }
  } // End of loop on vertices



  /*
   *  Prepare face triangulation
   */
  /* Count max nb of vertices per face */
  PDM_l_num_t n_vtx_face, n_vtx_face_max = 0;
  for (int iface = 0; iface < n_face; iface++) {
    n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];
    if (n_vtx_face > n_vtx_face_max) {
      n_vtx_face_max = n_vtx_face;
    }
  }

  /*
   *  Loop on faces
   */
  PDM_l_num_t n_tri;
  PDM_l_num_t _tri_vtx[3];
  PDM_l_num_t *tri_vtx = malloc (sizeof(PDM_l_num_t) * (n_vtx_face_max - 2)*3);
  PDM_triangulate_state_t *state = PDM_triangulate_state_create (n_vtx_face_max);

  for (int iface = 0; iface < n_face; iface++) {

    const PDM_l_num_t *_face_vtx = face_vtx + face_vtx_idx[iface];

    /*
     *  Triangulate current face
     */
    n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

    /* Triangular face */
    if (n_vtx_face == 3) {
      n_tri = 1;
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        tri_vtx[ivtx] = _face_vtx[ivtx];
      }
    }

    /* Quadrilateral face */
    else if (n_vtx_face == 4) {
      n_tri = PDM_triangulate_quadrangle (3,
                                          vtx_coord,
                                          NULL,
                                          _face_vtx,
                                          tri_vtx);
    }

    /* Polygonal face */
    else {
      n_tri = PDM_triangulate_polygon(3,
                                      n_vtx_face,
                                      vtx_coord,
                                      NULL,
                                      _face_vtx,
                                      PDM_TRIANGULATE_MESH_DEF,
                                      tri_vtx,
                                      state);
    }


    /* Loop on triangles */
    for (int itri = 0; itri < n_tri; itri++) {
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        _tri_vtx[ivtx] = tri_vtx[3*itri + ivtx] - 1;
      }

      if (face_orientation[iface] < 0) {
        PDM_l_num_t tmp = _tri_vtx[0];
        _tri_vtx[0] = _tri_vtx[2];
        _tri_vtx[2] = tmp;
      }

      double l[3] = {0., 0., 0.}, theta[3], sint[3], h = 0.;
      for (int i = 0; i < 3; i++) {
        int ip = (i+1)%3;
        int im = (i+2)%3;
        for (int idim = 0; idim < 3; idim++) {
          double delta = u[3*_tri_vtx[ip] + idim] - u[3*_tri_vtx[im] + idim];
          l[i] += delta * delta;
        }
        l[i] = sqrt(l[i]);

        theta[i] = asin(0.5 * l[i]);
        h += theta[i];
        theta[i] *= 2.;
        sint[i] = sin(theta[i]);
      }

      if (M_PI - h < eps) {
        /*
         *  point lies on current tirangle, use 2D barycentric coordinates
         */

        /* Triangular face */
        if (n_tri == 1 || LOCATE_ON_TRIANGLES) {
          for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
            mean_value_coord[ivtx] = 0.;
          }

          double sum_w = 0.;
          for (int i = 0; i < 3; i++) {
            int ivtx = _tri_vtx[i];
            int ip = _tri_vtx[(i+1)%3];
            int im = _tri_vtx[(i+2)%3];
            double w = sint[i] * d[ip] * d[im];
            sum_w += w;
            mean_value_coord[ivtx] = w;
          }

          for (int i = 0; i < 3; i++) {
            mean_value_coord[_tri_vtx[i]] /= sum_w;
          }
        }

        /* Polygonal face */
        else {
          n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];
          double *face_coord = malloc (sizeof(double) * n_vtx_face * 3);
          double *mean_value_coord_face = malloc (sizeof(double) * n_vtx_face);

          for (int j = 0; j < n_vtx_face; j++) {
            int id_vtx = face_vtx[face_vtx_idx[iface] + j] - 1;
            for (int idim = 0; idim < 3; idim++) {
              face_coord[3*j + idim] = vtx_coord[3*id_vtx + idim];
            }
          }

          PDM_mean_value_coordinates_polygon_3d (n_vtx_face,
                                                 face_coord,
                                                 1,
                                                 pt_coord,
                                                 mean_value_coord_face);

          for (int j = 0; j < n_vtx; j++) {
            mean_value_coord[j] = 0.;
          }

          for (int j = 0; j < n_vtx_face; j++) {
            int id_vtx = face_vtx[face_vtx_idx[iface] + j] - 1;
            mean_value_coord[id_vtx] = mean_value_coord_face[j];
          }

          free (face_coord);
          free (mean_value_coord_face);
        }

        free (u);
        free (d);
        free (tri_vtx);

        if (DEBUG) {
          printf("point located on face %d (M_PI - h = %g)\n", iface, M_PI - h);
        }
        return;
      }

      double c[3], s[3];
      double det = _determinant_3x3 ((u + 3*_tri_vtx[0]),
                                     (u + 3*_tri_vtx[1]),
                                     (u + 3*_tri_vtx[2]));
      PDM_bool_t ignore_triangle = PDM_FALSE;
      for (int i = 0; i < 3; i++) {
        c[i] = 2. * sin(h) * sin(h - theta[i]) / (sint[(i+1)%3] * sint[(i+2)%3]) - 1.;
        s[i] = sqrt(1. - c[i]*c[i]);

        if (s[i] < eps) {
          /* point lies outside current triangle on the same plane, ignore current triangle */
          ignore_triangle = PDM_TRUE;
          break;
        }

        if (det < 0.) {
          s[i] = -s[i];
        }
      }

      if (ignore_triangle == PDM_TRUE) {
        if (DEBUG) {
          printf("ignore triangle %d of face %d\n", itri, iface);
        }
        continue;
      }

      for (int i = 0; i < 3; i++) {
        int ip = (i+1)%3;
        int im = (i+2)%3;
        double w = (theta[i] - c[ip]*theta[im] - c[im]*theta[ip]) / (sint[ip] * s[im]);

        mean_value_coord[_tri_vtx[i]] += w;
      }

    } // End of loop on triangles

  } // End of loop on faces


  /* Normalize */
  double sum = 0.;
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    mean_value_coord[ivtx] /= d[ivtx];
    sum += mean_value_coord[ivtx];
  }

  if (fabs(sum) > 1.e-15) {
    sum = 1. / sum;
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      mean_value_coord[ivtx] *= sum;
    }
  }

  free (u);
  free (d);
  free (tri_vtx);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
