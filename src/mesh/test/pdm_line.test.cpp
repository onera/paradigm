#include "doctest/extensions/doctest_mpi.h"

#include <limits.h>
#include <float.h>
#include <math.h>
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_line.h"
#include "pdm_logging.h"


static const double tol = 1e-14;

TEST_CASE("PDM_ray_segment_intersection_2d - Valid intersection") {
    double u, v;
    double intersection_coord[3];
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {0.0, 1.0};
    const double b2[2] = {1.0, 0.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v, intersection_coord) == PDM_LINE_INTERSECT_YES);
    CHECK(u >= 0);  // Vérifie que le point d'intersection est sur la demi-droite
    CHECK(v >= 0);  // Vérifie que le point d'intersection est sur le segment
    CHECK(v <= 1);  // Vérifie que le point d'intersection est à l'intérieur du segment
    CHECK(fabs(intersection_coord[0] - 0.5) < tol);
    CHECK(fabs(intersection_coord[1] - 0.5) < tol);
}

TEST_CASE("PDM_ray_segment_intersection_2d - Parallel not collinear") {
    double u, v;
    double intersection_coord[3];
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {0.0, -1.0};
    const double b2[2] = {1.0, 0.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v, intersection_coord) == PDM_LINE_INTERSECT_NO);
}

TEST_CASE("PDM_ray_segment_intersection_2d - Intersection outside ray") {
    double u, v;
    double intersection_coord[3];
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {-1.0, 0.0};
    const double b2[2] = {-1.0, -2.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v, intersection_coord) == PDM_LINE_INTERSECT_NO);
}



TEST_CASE("PDM_ray_segment_intersection_2d - Intersection outside segment") {
    double u, v;
    double intersection_coord[3];
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {1.0, 0.0};
    const double b2[2] = {2.0, -1.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v, intersection_coord) == PDM_LINE_INTERSECT_NO);
}

TEST_CASE("PDM_ray_segment_intersection_2d - Colinear with intersection") {
    double u, v;
    double intersection_coord[3];
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {2.0, 2.0};
    const double b1[2] = {1.0, 1.0};
    const double b2[2] = {3.0, 3.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v, intersection_coord) == PDM_LINE_INTERSECT_ON_LINE);
    CHECK(fabs(u - 0.5) < tol); // Vérifie que le point d'intersection est sur la demi-droite
    CHECK(fabs(v)       < tol); // Vérifie que le point d'intersection est sur le point b1
    CHECK(fabs(intersection_coord[0] - b1[0]) < tol);
    CHECK(fabs(intersection_coord[1] - b1[1]) < tol);
}

TEST_CASE("PDM_ray_segment_intersection_2d - Colinear without intersection") { // Change name
    double u, v;
    double intersection_coord[3];
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {2.0, 2.0};
    const double b2[2] = {3.0, 3.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v, intersection_coord) == PDM_LINE_INTERSECT_ON_LINE);
    CHECK(fabs(u - 2) < tol); // Vérifie que le point d'intersection est sur la demi-droite
    CHECK(fabs(v)     < tol); // Vérifie que le point d'intersection est sur le point b1
    CHECK(fabs(intersection_coord[0] - b1[0]) < tol);
    CHECK(fabs(intersection_coord[1] - b1[1]) < tol);
}


TEST_CASE("PDM_ray_segment_intersection_2d - Intersection outside ray or segment") { // Change name
    double u, v;
    double intersection_coord[3];
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {0.5, 0.5};
    const double b2[2] = {-0.5, -0.5};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v, intersection_coord) == PDM_LINE_INTERSECT_ON_LINE);
    CHECK(fabs(u)       < tol); // Vérifie que le point d'intersection est sur le point a1
    CHECK(fabs(v - 0.5) < tol); // Vérifie que le point d'intersection est sur le segment
    CHECK(fabs(intersection_coord[0] - a1[0]) < tol);
    CHECK(fabs(intersection_coord[1] - a1[1]) < tol);
}

TEST_CASE("PDM_ray_segment_intersection_2d - Ray and segment coincident") {
    double u, v;
    double intersection_coord[3];
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {0.0, 0.0};
    const double b2[2] = {1.0, 1.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v, intersection_coord) == PDM_LINE_INTERSECT_ON_LINE);
    CHECK(fabs(u) < tol); // Vérifie que le point d'intersection est sur le point a1
    CHECK(fabs(v) < tol); // Vérifie que le point d'intersection est sur le point b1
    CHECK(fabs(intersection_coord[0] - a1[0]) < tol);
    CHECK(fabs(intersection_coord[1] - a1[1]) < tol);
}

TEST_CASE("PDM_ray_segment_intersection_2d - Segment is a point") {
    double u, v;
    double intersection_coord[3];
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {0.5, 0.5};
    const double b2[2] = {0.5, 0.5};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v, intersection_coord) == PDM_LINE_INTERSECT_ON_LINE);
    CHECK(fabs(u - 0.5) < tol); // Vérifie que le point d'intersection est sur le point a1
    CHECK(fabs(v) < tol); // Vérifie que le point d'intersection est sur le point b1
    CHECK(fabs(intersection_coord[0] - 0.5) < tol);
    CHECK(fabs(intersection_coord[1] - 0.5) < tol);
}

TEST_CASE("PDM_ray_segment_intersection_2d - Ray is a point") {
    double u, v;
    double intersection_coord[3];
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {0.0, 0.0};
    const double b1[2] = {0.5, 0.5};
    const double b2[2] = {-0.5, 0.5};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v, intersection_coord) == PDM_LINE_INTERSECT_UNDEF);
}


TEST_CASE("PDM_ray_segment_intersection_2d - Segment touches ray") {
    double u, v;
    double intersection_coord[3];
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {0.0, 1.0};
    const double b2[2] = {1.0, 1.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v, intersection_coord) == PDM_LINE_INTERSECT_YES);
    CHECK(fabs(u - 1) < tol); // Vérifie que le point d'intersection est sur le point a2
    CHECK(fabs(v - 1) < tol); // Vérifie que le point d'intersection est sur le point b2
    CHECK(fabs(intersection_coord[0] - a2[0]) < tol);
    CHECK(fabs(intersection_coord[1] - a2[1]) < tol);
}
