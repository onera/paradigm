#include "doctest/extensions/doctest_mpi.h"

#include <limits.h>
#include <float.h>
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_line.h"
#include "pdm_logging.h"


TEST_CASE("PDM_ray_segment_intersection_2d - Valid intersection") {
    double u, v;
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {0.0, 1.0};
    const double b2[2] = {1.0, 0.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v) == PDM_LINE_INTERSECT_YES);
    CHECK(u >= 0);  // Vérifie que le point d'intersection est sur la demi-droite
    CHECK(v >= 0);  // Vérifie que le point d'intersection est sur le segment
    CHECK(v <= 1);  // Vérifie que le point d'intersection est à l'intérieur du segment
}

TEST_CASE("PDM_ray_segment_intersection_2d - No intersection") {
    double u, v;
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {2.0, 2.0};  // Segment parallèle mais décalé
    const double b2[2] = {3.0, 3.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v) == PDM_LINE_INTERSECT_NO);
}

TEST_CASE("PDM_ray_segment_intersection_2d - Intersection outside segment") {
    double u, v;
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {1.0, 0.0};
    const double b2[2] = {2.0, -1.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v) == PDM_LINE_INTERSECT_NO);
}

TEST_CASE("PDM_ray_segment_intersection_2d - Colinear with intersection") {
    double u, v;
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {2.0, 2.0};
    const double b1[2] = {1.0, 1.0};
    const double b2[2] = {3.0, 3.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v) == PDM_LINE_INTERSECT_ON_LINE);
}

TEST_CASE("PDM_ray_segment_intersection_2d - Colinear without intersection") {
    double u, v;
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {2.0, 2.0};
    const double b2[2] = {3.0, 3.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v) == PDM_LINE_INTERSECT_NO);
}


TEST_CASE("PDM_ray_segment_intersection_2d - Intersection outside ray or segment") {
    double u, v;
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {0.5, 0.5};
    const double b2[2] = {-0.5, -0.5};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v) == PDM_LINE_INTERSECT_NO);
}

TEST_CASE("PDM_ray_segment_intersection_2d - Ray and segment coincident") {
    double u, v;
    const double a1[2] = {0.0, 0.0};
    const double a2[2] = {1.0, 1.0};
    const double b1[2] = {0.0, 0.0};
    const double b2[2] = {1.0, 1.0};

    CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v) == PDM_LINE_INTERSECT_ON_LINE);
}

// TEST_CASE("PDM_ray_segment_intersection_2d - Segment is a point") {
//     double u, v;
//     const double a1[2] = {0.0, 0.0};
//     const double a2[2] = {1.0, 1.0};
//     const double b1[2] = {0.5, 0.5};
//     const double b2[2] = {0.5, 0.5};

//     CHECK(PDM_ray_segment_intersection_2d(a1, a2, b1, b2, &u, &v) == PDM_LINE_INTERSECT_YES);
// }
