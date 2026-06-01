#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_quaternion.h"
#include "pdm_quaternion_priv.h"

static const double EPS = 8*__DBL_EPSILON__;
static const double DEG2RAD = M_PI/180.;
// static const double RAD2DEG = 180./M_PI;

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_set", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(0.,1.,2.,3.,&q);
    CHECK(PDM_quaternion_equal(&q,0.,1.,2.,3.,EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_equal", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(0.,1.,2.,3.,&q);
    CHECK( PDM_quaternion_equal(&q,0.,1.,2.,3.,EPS));
    CHECK(!PDM_quaternion_equal(&q,1.,1.,2.,3.,EPS));
    CHECK(!PDM_quaternion_equal(&q,0.,0.,2.,3.,EPS));
    CHECK(!PDM_quaternion_equal(&q,0.,1.,0.,3.,EPS));
    CHECK(!PDM_quaternion_equal(&q,0.,1.,2.,0.,EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_set_identity", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(0.,1.,2.,3.,&q);
    PDM_quaternion_set_identity(&q);
    CHECK(PDM_quaternion_equal(&q,1.,0.,0.,0.,EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_conjugate", 1) {
    PDM_quaternion q;
    PDM_quaternion qt;
    PDM_quaternion_set(0.,1.,2.,3.,&q);
    PDM_quaternion_set(3.,2.,1.,0.,&qt);
    PDM_quaternion_conjugate(&q,&qt);
    CHECK(PDM_quaternion_equal(&qt,q.q[0],-q.q[1],-q.q[2],-q.q[3],EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_norm", 1) {
    PDM_quaternion q;
    double norm;
    PDM_quaternion_set(3.,4.,0.,0.,&q);
    norm = PDM_quaternion_norm(&q);
    CHECK(PDM_ABS(norm-5.)<__DBL_EPSILON__);
    PDM_quaternion_set(0.,3.,4.,0.,&q);
    norm = PDM_quaternion_norm(&q);
    CHECK(PDM_ABS(norm-5.)<__DBL_EPSILON__);
    PDM_quaternion_set(0.,0.,3.,4.,&q);
    norm = PDM_quaternion_norm(&q);
    CHECK(PDM_ABS(norm-5.)<__DBL_EPSILON__);
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_normalize", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(0.,3.,4.,0.,&q);
    PDM_quaternion_normalize(&q);
    double norm = PDM_quaternion_norm(&q);
    CHECK(PDM_ABS(norm-1.)<__DBL_EPSILON__);

    PDM_quaternion_set(2.,0.,0.,0.,&q);
    PDM_quaternion_normalize(&q);
    // PDM_quaternion_print(&q);
    CHECK(PDM_quaternion_equal(&q,1.,0.,0.,0.,EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_compose", 1) {
    PDM_quaternion q1;
    PDM_quaternion q2;
    PDM_quaternion q_out;
    PDM_quaternion_set(0.,0.,0.,0.,&q1);
    PDM_quaternion_set(1.,2.,3.,4.,&q2);
    PDM_quaternion_set(0.,0.,0.,0.,&q_out);
    PDM_quaternion_set_identity(&q1);
    PDM_quaternion_compose(&q1,&q2,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&q2,EPS));

    PDM_quaternion_compose(&q2,&q1,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&q2,EPS));

}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_compose 2", 1) {
    PDM_quaternion q1;
    PDM_quaternion qi;
    PDM_quaternion qj;
    PDM_quaternion qk;
    PDM_quaternion q_out;
    PDM_quaternion_set(1.,0.,0.,0.,&q1);
    PDM_quaternion_set(0.,1.,0.,0.,&qi);
    PDM_quaternion_set(0.,0.,1.,0.,&qj);
    PDM_quaternion_set(0.,0.,0.,1.,&qk);
    PDM_quaternion_set(0.,0.,0.,0.,&q_out);

    PDM_quaternion_compose(&q1,&q1,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&q1,EPS));
    PDM_quaternion_compose(&q1,&qi,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&qi,EPS));
    PDM_quaternion_compose(&q1,&qj,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&qj,EPS));
    PDM_quaternion_compose(&q1,&qk,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&qk,EPS));

    PDM_quaternion_compose(&qi,&q1,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&qi,EPS));
    PDM_quaternion_compose(&qi,&qi,&q_out);
    CHECK(PDM_quaternion_equal(&q_out,-q1.q[0],-q1.q[1],-q1.q[2],-q1.q[3],EPS));
    PDM_quaternion_compose(&qi,&qj,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&qk,EPS));
    PDM_quaternion_compose(&qi,&qk,&q_out);
    CHECK(PDM_quaternion_equal(&q_out,-qj.q[0],-qj.q[1],-qj.q[2],-qj.q[3],EPS));

    PDM_quaternion_compose(&qj,&q1,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&qj,EPS));
    PDM_quaternion_compose(&qj,&qi,&q_out);
    CHECK(PDM_quaternion_equal(&q_out,-qk.q[0],-qk.q[1],-qk.q[2],-qk.q[3],EPS));
    PDM_quaternion_compose(&qj,&qj,&q_out);
    CHECK(PDM_quaternion_equal(&q_out,-q1.q[0],-q1.q[1],-q1.q[2],-q1.q[3],EPS));
    PDM_quaternion_compose(&qj,&qk,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&qi,EPS));

    PDM_quaternion_compose(&qk,&q1,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&qk,EPS));
    PDM_quaternion_compose(&qk,&qi,&q_out);
    CHECK(PDM_quaternion_equal_quaternion(&q_out,&qj,EPS));
    PDM_quaternion_compose(&qk,&qj,&q_out);
    CHECK(PDM_quaternion_equal(&q_out,-qi.q[0],-qi.q[1],-qi.q[2],-qi.q[3],EPS));
    PDM_quaternion_compose(&qk,&qk,&q_out);
    CHECK(PDM_quaternion_equal(&q_out,-q1.q[0],-q1.q[1],-q1.q[2],-q1.q[3],EPS));

}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_rotate", 1) {
    // rotation of 120deg of axis 1,1,1 (i->k,k->j,j->i)
    PDM_quaternion q;
    PDM_quaternion_set(-.5,.5,.5,.5,&q);
    double vector[9] = {1.,0.,0.,
        0.,1.,0.,
        0.,0.,1.,
    };
    double out_vect[9];
    double exp_vect[9] = {0.,0.,1.,
        1.,0.,0.,
        0.,1.,0.};
        PDM_quaternion_rotate(&q,vector,3,out_vect);
        CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,9,EPS);

    PDM_quaternion qt;
    PDM_quaternion_set(-.5,-.5,-.5,-.5,&qt);
    double exp_vect_2[9] = {0.,1.,0.,
                            0.,0.,1.,
                            1.,0.,0.};
    PDM_quaternion_rotate(&qt,vector,3,out_vect);
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect_2,9,EPS);
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_rotate 2", 1) {

    // For checks: https://www.andre-gaschler.com/rotationconverter/
    PDM_quaternion q;
    double sin_45 = .5*sqrt(2.);
    double vector[3] = {1.,2.,3.};
    double out_vect[3];
    double exp_vect[3];

    // rotation of 90 deg around x-axis
    PDM_quaternion_set(sin_45,sin_45,0.,0.,&q);
    exp_vect[0] =  1.;
    exp_vect[1] = -3.;
    exp_vect[2] =  2.;
    PDM_quaternion_rotate(&q,vector,1,out_vect);
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);

    // rotation of 90 deg around y-axis
    PDM_quaternion_set(sin_45,0.,sin_45,0.,&q);
    exp_vect[0] =  3.;
    exp_vect[1] =  2.;
    exp_vect[2] = -1.;
    PDM_quaternion_rotate(&q,vector,1,out_vect);
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);

    // rotation of 90 deg around z-axis
    PDM_quaternion_set(sin_45,0.,0.,sin_45,&q);
    exp_vect[0] = -2.;
    exp_vect[1] =  1.;
    exp_vect[2] =  3.;
    PDM_quaternion_rotate(&q,vector,1,out_vect);
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);

    // arbitrary rotation of 37 deg around a (1,2,3)-axis
    double cos_18_5 = cos(37.*0.5*DEG2RAD);
    double sin_18_5 = sin(37.*0.5*DEG2RAD);
    double axis[3] = {1.,2.,3.};
    double ax_inv_norm = 1./sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
    axis[0] *= ax_inv_norm;
    axis[1] *= ax_inv_norm;
    axis[2] *= ax_inv_norm;
    PDM_quaternion_set(cos_18_5,axis[0]*sin_18_5,axis[1]*sin_18_5,axis[2]*sin_18_5,&q);
    vector[0] = 1.;
    vector[1] = 0.;
    vector[2] = 0.;
    exp_vect[0] =  0.8130186879010576;
    exp_vect[1] =  0.5112918471750423;
    exp_vect[2] = -0.2785341274170473;
    PDM_quaternion_rotate(&q,vector,1,out_vect);
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);
    vector[0] = 0.;
    vector[1] = 1.;
    vector[2] = 0.;
    exp_vect[0] = -0.45375913575998306;
    exp_vect[1] =  0.856168221462352;
    exp_vect[2] =  0.24714089761175967;
    PDM_quaternion_rotate(&q,vector,1,out_vect);
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);
    vector[0] = 0.;
    vector[1] = 0.;
    vector[2] = 1.;
    exp_vect[0] =  0.36483319453963614;
    exp_vect[1] = -0.07454276336658205;
    exp_vect[2] =  0.928084110731176;
    PDM_quaternion_rotate(&q,vector,1,out_vect);
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_two_vectors", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(0.,0.,0.,0.,&q);
    double v1[3] = {1.,0.,0.};
    double v2[3] = {0.,1.,0.};
    PDM_quaternion_from_two_vectors(v1,v2,&q);
    double sin_45 = .5*sqrt(2.);
    CHECK(PDM_quaternion_equal(&q,sin_45,0.,0.,sin_45,EPS));

    double v3[3];
    PDM_quaternion_rotate(&q,v1,1,v3);
    CHECK_EQ_C_ARRAY_FLOAT(v3,v2,3,EPS);

    // arbitrary test case
    v1[0] =  1.;
    v1[1] = -2.;
    v1[2] =  3.;
    v2[0] = -4.;
    v2[1] =  5.;
    v2[2] = -6.;
    PDM_quaternion_from_two_vectors(v1,v2,&q);
    double v1_inv_norm = 1./sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
    double v2_inv_norm = 1./sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
    v1[0] *= v1_inv_norm;
    v1[1] *= v1_inv_norm;
    v1[2] *= v1_inv_norm;
    v2[0] *= v2_inv_norm;
    v2[1] *= v2_inv_norm;
    v2[2] *= v2_inv_norm;
    PDM_quaternion_rotate(&q,v1,1,v3);
    CHECK_EQ_C_ARRAY_FLOAT(v3,v2,3,EPS);

    // parallel vectors
    v1[0] =  1.;
    v1[1] = -2.;
    v1[2] =  3.;
    v2[0] =  2.;
    v2[1] = -4.;
    v2[2] =  6.;
    PDM_quaternion_from_two_vectors(v1,v2,&q);
    CHECK(PDM_quaternion_equal(&q,1.,0.,0.,0.,EPS)); // expects identity
    v1_inv_norm = 1./sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
    v2_inv_norm = 1./sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
    v1[0] *= v1_inv_norm;
    v1[1] *= v1_inv_norm;
    v1[2] *= v1_inv_norm;
    v2[0] *= v2_inv_norm;
    v2[1] *= v2_inv_norm;
    v2[2] *= v2_inv_norm;
    PDM_quaternion_rotate(&q,v1,1,v3);
    CHECK_EQ_C_ARRAY_FLOAT(v3,v2,3,EPS);

    // opposite vectors v1 is z axis
    v1[0] =  0.;
    v1[1] =  0.;
    v1[2] = -3.;
    v2[0] =  0.;
    v2[1] =  0.;
    v2[2] =  1.;
    PDM_quaternion_from_two_vectors(v1,v2,&q);
    v1_inv_norm = 1./sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
    v2_inv_norm = 1./sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
    v1[0] *= v1_inv_norm;
    v1[1] *= v1_inv_norm;
    v1[2] *= v1_inv_norm;
    v2[0] *= v2_inv_norm;
    v2[1] *= v2_inv_norm;
    v2[2] *= v2_inv_norm;
    PDM_quaternion_rotate(&q,v1,1,v3);
    CHECK_EQ_C_ARRAY_FLOAT(v3,v2,3,EPS);

    // opposite vectors v1 is not z axis
    v1[0] =  1.;
    v1[1] =  2.;
    v1[2] = -3.;
    v2[0] = -2.;
    v2[1] = -4.;
    v2[2] =  6.;
    PDM_quaternion_from_two_vectors(v1,v2,&q);
    v1_inv_norm = 1./sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
    v2_inv_norm = 1./sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
    v1[0] *= v1_inv_norm;
    v1[1] *= v1_inv_norm;
    v1[2] *= v1_inv_norm;
    v2[0] *= v2_inv_norm;
    v2[1] *= v2_inv_norm;
    v2[2] *= v2_inv_norm;
    PDM_quaternion_rotate(&q,v1,1,v3);
    CHECK_EQ_C_ARRAY_FLOAT(v3,v2,3,EPS);
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_axis_aligned_rotation", 1) {
    PDM_quaternion q;
    double ang = 37*DEG2RAD;
    double cos_a = cos(ang*0.5);
    double sin_a = sin(ang*0.5);

    PDM_quaternion_from_axis_aligned_rotation(ang,0,&q);
    CHECK(PDM_quaternion_equal(&q,cos_a,sin_a,0.,0.,EPS));
    PDM_quaternion_from_axis_aligned_rotation(ang,1,&q);
    CHECK(PDM_quaternion_equal(&q,cos_a,0.,sin_a,0.,EPS));
    PDM_quaternion_from_axis_aligned_rotation(ang,2,&q);
    CHECK(PDM_quaternion_equal(&q,cos_a,0.,0.,sin_a,EPS));

}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_axis_aligned_symmetry", 1) {
    PDM_quaternion q;
    double cos_a = 0.;
    double sin_a = 1.;

    PDM_quaternion_from_axis_aligned_symmetry(0,&q);
    CHECK(PDM_quaternion_equal(&q,cos_a,sin_a,0.,0.,EPS));
    PDM_quaternion_from_axis_aligned_symmetry(1,&q);
    CHECK(PDM_quaternion_equal(&q,cos_a,0.,sin_a,0.,EPS));
    PDM_quaternion_from_axis_aligned_symmetry(2,&q);
    CHECK(PDM_quaternion_equal(&q,cos_a,0.,0.,sin_a,EPS));

}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_axis_angle", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(0.,0.,0.,0.,&q);
    double axis[3] = {1.,1.,1.};
    double angle = 120.*DEG2RAD;
    PDM_quaternion_from_axis_angle(axis,angle,&q);

    CHECK(PDM_quaternion_equal(&q,.5,.5,.5,.5,EPS));

    double axis2[3] = {-1.,-1.,-1.};
    angle = -120.*DEG2RAD;
    PDM_quaternion_from_axis_angle(axis2,angle,&q);

    CHECK(PDM_quaternion_equal(&q,.5,.5,.5,.5,EPS));

    angle = 37*DEG2RAD;
    axis[0] = -1.;
    axis[1] =  2.;
    axis[2] = -3.;
    PDM_quaternion_from_axis_angle(axis,angle,&q);
    CHECK(PDM_quaternion_equal(&q,9.4832365520619932e-01,-8.4803236535420032e-02, 1.6960647307084006e-01,-2.5440970960626014e-01,EPS));

}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_axis_angle 2", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(0.,0.,0.,0.,&q);
    double angle;
    double axes[75] = {
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 1.0, 1.0,
        1.0, -1.0, 1.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 1.0, 1.0,
        1.0, -1.0, 1.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 1.0, 1.0,
        1.0, -1.0, 1.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 1.0, 1.0,
        1.0, -1.0, 1.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 1.0, 1.0,
        1.0, -1.0, 1.0
    };
    double w[25] = {
        0.8660254037844387,
        0.8660254037844387,
        0.8660254037844387,
        0.8660254037844386,
        0.8660254037844386,
        0.9659258262890683,
        0.9659258262890683,
        0.9659258262890683,
        0.9659258262890683,
        0.9659258262890683,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        0.9659258262890683,
        0.9659258262890683,
        0.9659258262890683,
        0.9659258262890683,
        0.9659258262890683,
        0.8660254037844387,
        0.8660254037844387,
        0.8660254037844387,
        0.8660254037844386,
        0.8660254037844386
    };
    double v0[25] = {
       -0.49999999999999994,
       -0.0,
       -0.0,
       -0.28867513459481287,
       -0.28867513459481287,
       -0.25881904510252074,
       -0.0,
       -0.0,
       -0.14942924536134228,
       -0.14942924536134228,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.25881904510252074,
        0.0,
        0.0,
        0.14942924536134228,
        0.14942924536134228,
        0.49999999999999994,
        0.0,
        0.0,
        0.28867513459481287,
        0.28867513459481287
    };
    double v1[25] = {
       -0.0,
       -0.49999999999999994,
       -0.0,
       -0.28867513459481287,
        0.28867513459481287,
       -0.0,
       -0.25881904510252074,
       -0.0,
       -0.14942924536134228,
        0.14942924536134228,
        0.0,
        0.0,
        0.0,
        0.0,
       -0.0,
        0.0,
        0.25881904510252074,
        0.0,
        0.14942924536134228,
       -0.14942924536134228,
        0.0,
        0.49999999999999994,
        0.0,
        0.28867513459481287,
       -0.28867513459481287
    };
    double v2[25] = {
       -0.0,
       -0.0,
       -0.49999999999999994,
       -0.28867513459481287,
       -0.28867513459481287,
       -0.0,
       -0.0,
       -0.25881904510252074,
       -0.14942924536134228,
       -0.14942924536134228,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.25881904510252074,
        0.14942924536134228,
        0.14942924536134228,
        0.0,
        0.0,
        0.49999999999999994,
        0.28867513459481287,
        0.28867513459481287
    };

    for (int i = 0; i < 5; i++) {
        angle = (-60.+30.*i)*DEG2RAD;
        for (int j = 0; j < 5; j++) {
            PDM_quaternion_from_axis_angle(&axes[3*(5*i+j)],angle,&q);
            CHECK(PDM_quaternion_equal(&q,w[5*i+j],v0[5*i+j],v1[5*i+j],v2[5*i+j],EPS));
        }
    }
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_euler_angles intrinsic [2,1,0]", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(0.,0.,0.,0.,&q);
    double sin_45 = .5*sqrt(2.);
    PDM_bool_t intrinsic = PDM_TRUE;
    int order[3] = {2,1,0};

    PDM_quaternion_from_euler_angles(.5*M_PI,0.,0.,order,intrinsic,&q);
    CHECK(PDM_quaternion_equal(&q,sin_45,sin_45,0.,0.,EPS));
    PDM_quaternion_from_euler_angles(0.,0.5*M_PI,0.,order,intrinsic,&q);
    CHECK(PDM_quaternion_equal(&q,sin_45,0.,sin_45,0.,EPS));
    PDM_quaternion_from_euler_angles(0.,0.,0.5*M_PI,order,intrinsic,&q);
    CHECK(PDM_quaternion_equal(&q,sin_45,0.,0.,sin_45,EPS));
    PDM_quaternion_from_euler_angles(0.5*M_PI,0.5*M_PI,0.,order,intrinsic,&q);
    CHECK(PDM_quaternion_equal(&q,.5,.5,.5,-.5,EPS));
    PDM_quaternion_from_euler_angles(0.5*M_PI,0.,0.5*M_PI,order,intrinsic,&q);
    CHECK(PDM_quaternion_equal(&q,.5,.5,.5,.5,EPS));
    PDM_quaternion_from_euler_angles(0.,0.5*M_PI,0.5*M_PI,order,intrinsic,&q);
    CHECK(PDM_quaternion_equal(&q,.5,-.5,.5,.5,EPS));

}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_euler_angles 2 intrinsic [2,1,0]", 1) {
    PDM_quaternion _q;
    PDM_quaternion* q = &_q;
    PDM_quaternion_set(0.,0.,0.,0.,q);
    PDM_bool_t intrinsic = PDM_TRUE;
    int order[3] = {2,1,0};
    double ang_x,ang_y,ang_z;
    double w[27] = {
        0.9879654343559628,
        0.9924038765061041,
        0.9892895259261897,
        0.9924038765061041,
        0.9961946980917455,
        0.9924038765061041,
        0.9892895259261897,
        0.9924038765061041,
        0.9879654343559628,
        0.9924038765061041,
        0.9961946980917455,
        0.9924038765061041,
        0.9961946980917455,
        1.0,
        0.9961946980917455,
        0.9924038765061041,
        0.9961946980917455,
        0.9924038765061041,
        0.9892895259261897,
        0.9924038765061041,
        0.9879654343559628,
        0.9924038765061041,
        0.9961946980917455,
        0.9924038765061041,
        0.9879654343559628,
        0.9924038765061041,
        0.9892895259261897
    };
    double v0[27] = {
       -0.09406091491321404,
       -0.007596123493895969,
        0.07892647901187543,
       -0.08682408883346517,
        0.0,
        0.08682408883346517,
       -0.07892647901187543,
        0.007596123493895969,
        0.09406091491321404,
       -0.08682408883346517,
        0.0,
        0.08682408883346517,
       -0.08715574274765817,
        0.0,
        0.08715574274765817,
       -0.08682408883346517,
        0.0,
        0.08682408883346517,
       -0.07892647901187543,
        0.007596123493895969,
        0.09406091491321404,
       -0.08682408883346517,
        0.0,
        0.08682408883346517,
       -0.09406091491321404,
       -0.007596123493895969,
        0.07892647901187543
    };
    double v1[27] = {
       -0.07892647901187541,
       -0.08682408883346517,
       -0.09406091491321403,
        0.007596123493895969,
        0.0,
       -0.007596123493895969,
        0.09406091491321403,
        0.08682408883346517,
        0.07892647901187541,
       -0.08682408883346517,
       -0.08715574274765817,
       -0.08682408883346517,
        0.0,
        0.0,
        0.0,
        0.08682408883346517,
        0.08715574274765817,
        0.08682408883346517,
       -0.09406091491321403,
       -0.08682408883346517,
       -0.07892647901187541,
       -0.007596123493895969,
        0.0,
        0.007596123493895969,
        0.07892647901187541,
        0.08682408883346517,
        0.09406091491321403
    };
    double v2[27] = {
       -0.09406091491321403,
       -0.08682408883346517,
       -0.07892647901187541,
       -0.08682408883346517,
       -0.08715574274765817,
       -0.08682408883346517,
       -0.07892647901187541,
       -0.08682408883346517,
       -0.09406091491321403,
       -0.007596123493895969,
        0.0,
        0.007596123493895969,
        0.0,
        0.0,
        0.0,
        0.007596123493895969,
        0.0,
       -0.007596123493895969,
        0.07892647901187541,
        0.08682408883346517,
        0.09406091491321403,
        0.08682408883346517,
        0.08715574274765817,
        0.08682408883346517,
        0.09406091491321403,
        0.08682408883346517,
        0.07892647901187541
    };
    for (int k = 0; k<3; k++) {
        ang_z = -10.+k*10.;
        for (int j = 0; j<3; j++) {
            ang_y = -10.+j*10.;
            for (int i = 0; i<3; i++) {
                ang_x = -10.+i*10.;
                PDM_quaternion_from_euler_angles(ang_x*DEG2RAD,ang_y*DEG2RAD,ang_z*DEG2RAD,order,intrinsic,q);
                CHECK(PDM_quaternion_equal(q,w[9*k+3*j+i],v0[9*k+3*j+i],v1[9*k+3*j+i],v2[9*k+3*j+i],EPS));
            }
        }
    }
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_euler_angles intrinsic/extrinsic ", 1) {
    PDM_quaternion q;
    PDM_quaternion q_expec;
    PDM_bool_t intrinsic = PDM_TRUE;
    PDM_bool_t extrinsic = PDM_FALSE;
    int order[3],rev_order[3];
    double ang_x = 15.*DEG2RAD;
    double ang_y = 25.*DEG2RAD;
    double ang_z = 35.*DEG2RAD;

    // XYZ
    order[0] = 0;
    order[1] = 1;
    order[2] = 2;
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&q);
    PDM_quaternion_set(
        9.1464902421112626e-01,1.8606208846158473e-01,1.6633655702976691e-01,3.1800976642617612e-01,
        &q_expec);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,rev_order,extrinsic,&q);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));

    // XZY
    order[0] = 0;
    order[1] = 2;
    order[2] = 1;
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&q);
    PDM_quaternion_set(
        9.3163952654103022e-01,5.7006410511950115e-02,1.6633655702976691e-01,3.1800976642617612e-01,
        &q_expec);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,rev_order,extrinsic,&q);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));

    // YXZ
    order[0] = 1;
    order[1] = 0;
    order[2] = 2;
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&q);
    PDM_quaternion_set(
        9.3163952654103022e-01,1.8606208846158473e-01,1.6633655702976691e-01,2.6412277754711250e-01,
        &q_expec);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,rev_order,extrinsic,&q);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));

    // YZX
    order[0] = 1;
    order[1] = 2;
    order[2] = 0;
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&q);
    PDM_quaternion_set(
        9.1464902421112626e-01,1.8606208846158473e-01,2.4297576037075486e-01,2.6412277754711250e-01,
        &q_expec);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,rev_order,extrinsic,&q);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));

    // ZXY
    order[0] = 2;
    order[1] = 0;
    order[2] = 1;
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&q);
    PDM_quaternion_set(
        9.1464902421112626e-01,5.7006410511950115e-02,2.4297576037075486e-01,3.1800976642617612e-01,
        &q_expec);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,rev_order,extrinsic,&q);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));

    // ZYX
    order[0] = 2;
    order[1] = 1;
    order[2] = 0;
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&q);
    PDM_quaternion_set(
        9.3163952654103022e-01,5.7006410511950115e-02,2.4297576037075486e-01,2.6412277754711250e-01,
        &q_expec);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,rev_order,extrinsic,&q);
    CHECK(PDM_quaternion_equal_quaternion(&q,&q_expec,EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_rotation_matrix", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(0.,0.,0.,0.,&q);
    double ang_x,ang_y,ang_z;
    PDM_bool_t intrinsic = PDM_TRUE;
    int order[3] = {2,1,0};
    double rot_mat[9];
    const double* r_mat = (const double*) rot_mat;
    double out_ang_x,out_ang_y,out_ang_z;

    for (int k = 0; k<3; k++) {
        ang_z = (-10.+k*10.)*DEG2RAD;
        for (int j = 0; j<3; j++) {
            ang_y = (-10.+j*10.)*DEG2RAD;
            for (int i = 0; i<3; i++) {
                ang_x = (-10.+i*10.)*DEG2RAD;
                rot_mat[3*0+0] = cos(ang_y)*cos(ang_z);
                rot_mat[3*1+0] = cos(ang_y)*sin(ang_z);
                rot_mat[3*2+0] = -sin(ang_y);
                rot_mat[3*0+1] = sin(ang_x)*sin(ang_y)*cos(ang_z)-cos(ang_x)*sin(ang_z);
                rot_mat[3*1+1] = sin(ang_x)*sin(ang_y)*sin(ang_z)+cos(ang_x)*cos(ang_z);
                rot_mat[3*2+1] = sin(ang_x)*cos(ang_y);
                rot_mat[3*0+2] = cos(ang_x)*sin(ang_y)*cos(ang_z)+sin(ang_x)*sin(ang_z);
                rot_mat[3*1+2] = cos(ang_x)*sin(ang_y)*sin(ang_z)-sin(ang_x)*cos(ang_z);
                rot_mat[3*2+2] = cos(ang_x)*cos(ang_y);

                PDM_quaternion_from_rotation_matrix(r_mat,&q);
                PDM_quaternion_to_euler_angles(&q,order,intrinsic,&out_ang_x,&out_ang_y,&out_ang_z);
                CHECK(PDM_ABS(ang_x-out_ang_x)<EPS);
                CHECK(PDM_ABS(ang_y-out_ang_y)<EPS);
                CHECK(PDM_ABS(ang_z-out_ang_z)<EPS);
            }
        }
    }
}
MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_rotation_matrix 2", 1) {
    PDM_quaternion q;

    // identity -> choice 3
    double rot_mat_0[9] = {
        1.,0.,0.,
        0.,1.,0.,
        0.,0.,1.
    };
    PDM_quaternion_from_rotation_matrix(rot_mat_0,&q);
    // printf("%s::%d\n",__FILE__,__LINE__);
    // PDM_quaternion_print(&q);
    CHECK(PDM_quaternion_equal(&q,1.,0.,0.,0.,EPS));

    // arbitrary -> choice 2
    double axis[3] = {1.,-2.,3.};
    double angle = 135*DEG2RAD;
    double rot_mat[9];
    PDM_quaternion_from_axis_angle(axis,angle,&q);
    // printf("%s::%d \n",__FILE__,__LINE__);
    // PDM_quaternion_print(&q);
    PDM_quaternion_to_rotation_matrix(&q,rot_mat);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[0*3+0],rot_mat[0*3+1],rot_mat[0*3+2]);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[1*3+0],rot_mat[1*3+1],rot_mat[1*3+2]);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[2*3+0],rot_mat[2*3+1],rot_mat[2*3+2]);
    PDM_quaternion_from_rotation_matrix(rot_mat,&q);
    // PDM_quaternion_print(&q);
    CHECK(PDM_quaternion_equal(&q,3.8268343236508995e-01,2.4691719123643649e-01,-4.9383438247287298e-01,7.4075157370930944e-01,EPS));

    // 120° around [1,1,1] => tr == 0 // choice 0
    double rot_mat_1[9] = {
        0.,0.,1.,
        1.,0.,0.,
        0.,1.,0.
    };
    // printf("%s::%d\n",__FILE__,__LINE__);
    PDM_quaternion_from_rotation_matrix(rot_mat_1,&q);
    // PDM_quaternion_print(&q);
    CHECK(PDM_quaternion_equal(&q,0.5,0.5,0.5,0.5,EPS));

    // 180° around [1,1,1] -> choice 3
    axis[0] = 1.;
    axis[1] = 1.;
    axis[2] = 1.;
    angle = 0.5*M_PI;
    PDM_quaternion_from_axis_angle(axis,angle,&q);
    // printf("%s::%d \n",__FILE__,__LINE__);
    // PDM_quaternion_print(&q);
    PDM_quaternion_to_rotation_matrix(&q,rot_mat);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[0*3+0],rot_mat[0*3+1],rot_mat[0*3+2]);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[1*3+0],rot_mat[1*3+1],rot_mat[1*3+2]);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[2*3+0],rot_mat[2*3+1],rot_mat[2*3+2]);
    PDM_quaternion_from_rotation_matrix(rot_mat,&q);
    // PDM_quaternion_print(&q);
    CHECK(PDM_quaternion_equal(&q,7.0710678118654746e-01,4.0824829046386296e-01,4.0824829046386296e-01,4.0824829046386296e-01,EPS));

    // 120° around [1,-2,1] -> choice 1
    axis[0] = 1.;
    axis[1] = -2.;
    axis[2] = 1.;
    angle = -120*DEG2RAD;
    PDM_quaternion_from_axis_angle(axis,angle,&q);
    // printf("%s::%d \n",__FILE__,__LINE__);
    PDM_quaternion_to_rotation_matrix(&q,rot_mat);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[0*3+0],rot_mat[0*3+1],rot_mat[0*3+2]);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[1*3+0],rot_mat[1*3+1],rot_mat[1*3+2]);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[2*3+0],rot_mat[2*3+1],rot_mat[2*3+2]);
    // PDM_quaternion_print(&q);
    PDM_quaternion_from_rotation_matrix(rot_mat,&q);
    CHECK(PDM_quaternion_equal(&q,5.0000000000000011e-01,-3.5355339059327379e-01,7.0710678118654757e-01,-3.5355339059327379e-01,EPS));

    // 120° around [1,1,-2] -> choice 2
    axis[0] = 1.;
    axis[1] = 1.;
    axis[2] = -2.;
    angle = -120*DEG2RAD;
    PDM_quaternion_from_axis_angle(axis,angle,&q);
    // printf("%s::%d \n",__FILE__,__LINE__);
    PDM_quaternion_to_rotation_matrix(&q,rot_mat);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[0*3+0],rot_mat[0*3+1],rot_mat[0*3+2]);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[1*3+0],rot_mat[1*3+1],rot_mat[1*3+2]);
    // printf("[%23.16e %23.16e %23.16e]\n",rot_mat[2*3+0],rot_mat[2*3+1],rot_mat[2*3+2]);
    PDM_quaternion_from_rotation_matrix(rot_mat,&q);
    // PDM_quaternion_print(&q);
    CHECK(PDM_quaternion_equal(&q,5.0000000000000011e-01,-3.5355339059327379e-01,-3.5355339059327379e-01,7.0710678118654735e-01,EPS));

}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_to_axis_angle", 1) {
    PDM_quaternion _q;
    PDM_quaternion* q = &_q;
    PDM_quaternion_set(0.,0.,0.,0.,q);
    double axis[3],expec_axis[3];
    double angle,out_angle;
    double inv_sqrt3 = 1./sqrt(3.);
    double axes[75] = {
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        inv_sqrt3, inv_sqrt3, inv_sqrt3,
        inv_sqrt3, -inv_sqrt3, inv_sqrt3,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        inv_sqrt3, inv_sqrt3, inv_sqrt3,
        inv_sqrt3, -inv_sqrt3, inv_sqrt3,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        inv_sqrt3, inv_sqrt3, inv_sqrt3,
        inv_sqrt3, -inv_sqrt3, inv_sqrt3,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        inv_sqrt3, inv_sqrt3, inv_sqrt3,
        inv_sqrt3, -inv_sqrt3, inv_sqrt3,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        inv_sqrt3, inv_sqrt3, inv_sqrt3,
        inv_sqrt3, -inv_sqrt3, inv_sqrt3
    };

    for (int k = 0; k < 5; k++) {
        angle = (-60.+30.*k)*DEG2RAD;
        for (int j = 0; j < 5; j++) {
            PDM_quaternion_from_axis_angle(&axes[3*(5*k+j)],angle,q);
            PDM_quaternion_to_axis_angle(q,axis,&out_angle);
            if (PDM_ABS(angle+out_angle)<EPS){ // results can be flipped (axis and angle)
                out_angle *= -1.;
                axis[0] *= -1.;
                axis[1] *= -1.;
                axis[2] *= -1.;
            }
            CHECK(PDM_ABS(angle-out_angle)<EPS);
            expec_axis[0] = axes[3*(5*k+j)+0];
            expec_axis[1] = axes[3*(5*k+j)+1];
            expec_axis[2] = axes[3*(5*k+j)+2];
            if (PDM_ABS(angle)>EPS){
                CHECK_EQ_C_ARRAY_FLOAT(expec_axis,axis,3,EPS);
            }
        }
    }
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_to_euler_angles", 1) {
    PDM_quaternion _q;
    PDM_quaternion* q = &_q;
    PDM_quaternion_set(0.,0.,0.,0.,q);
    PDM_bool_t intrinsic = PDM_TRUE;
    int order[3] = {2,1,0};
    double ang_x,ang_y,ang_z;
    double out_ang_x,out_ang_y,out_ang_z;
    for (int k = 0; k<3; k++) {
        ang_z = -10.+k*10.;
        for (int j = 0; j<3; j++) {
            ang_y = -10.+j*10.;
            for (int i = 0; i<3; i++) {
                ang_x = -10.+i*10.;
                PDM_quaternion_from_euler_angles(ang_x*DEG2RAD,ang_y*DEG2RAD,ang_z*DEG2RAD,order,intrinsic,q);
                PDM_quaternion_to_euler_angles(q,order,intrinsic,&out_ang_x,&out_ang_y,&out_ang_z);
                CHECK(PDM_ABS(ang_x*DEG2RAD-out_ang_x)<EPS);
                CHECK(PDM_ABS(ang_y*DEG2RAD-out_ang_y)<EPS);
                CHECK(PDM_ABS(ang_z*DEG2RAD-out_ang_z)<EPS);
            }
        }
    }
}
MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_to_euler_angles intrinsic/extrinsic", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(1.,-2.,3.,-4,&q);
    PDM_bool_t intrinsic = PDM_TRUE;
    PDM_bool_t extrinsic = PDM_FALSE;
    double ang_x,ang_y,ang_z;
    double exp_x,exp_y,exp_z;
    int order[3],rev_order[3];

    // XYZ
    order[0] = 0;
    order[1] = 1;
    order[2] = 2;
    PDM_quaternion_to_euler_angles(&q,order,intrinsic,&ang_x,&ang_y,&ang_z);
    exp_x = 1.3734007669450157e+00;
    exp_y = 8.2321197712587590e-01;
    exp_z = 2.9441970937399122e+00;
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    PDM_quaternion_to_euler_angles(&q,rev_order,extrinsic,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // XZY
    order[0] = 0;
    order[1] = 2;
    order[2] = 1;
    PDM_quaternion_to_euler_angles(&q,order,intrinsic,&ang_x,&ang_y,&ang_z);
    exp_x = -1.9138202672155999e+00;
    exp_y =  2.3086113869153615e+00;
    exp_z =  1.3373158940994179e-01;
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    PDM_quaternion_to_euler_angles(&q,rev_order,extrinsic,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // YXZ
    order[0] = 1;
    order[1] = 0;
    order[2] = 2;
    PDM_quaternion_to_euler_angles(&q,order,intrinsic,&ang_x,&ang_y,&ang_z);
    exp_x =  7.2972765622696656e-01;
    exp_y =  1.3909428270024184e+00;
    exp_z = -2.0344439357957027e+00;
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    PDM_quaternion_to_euler_angles(&q,rev_order,extrinsic,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // YZX
    order[0] = 1;
    order[1] = 2;
    order[2] = 0;
    PDM_quaternion_to_euler_angles(&q,order,intrinsic,&ang_x,&ang_y,&ang_z);
    exp_x =  2.0344439357957027e+00;
    exp_y = -2.6779450445889870e+00;
    exp_z = -7.2972765622696634e-01;
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    PDM_quaternion_to_euler_angles(&q,rev_order,extrinsic,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // ZXY
    order[0] = 2;
    order[1] = 0;
    order[2] = 1;
    PDM_quaternion_to_euler_angles(&q,order,intrinsic,&ang_x,&ang_y,&ang_z);
    exp_x = -1.2035883062370594e+00;
    exp_y = -1.1902899496825317e+00;
    exp_z =  2.7610862764774282e+00;
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    PDM_quaternion_to_euler_angles(&q,rev_order,extrinsic,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // ZYX
    order[0] = 2;
    order[1] = 1;
    order[2] = 0;
    PDM_quaternion_to_euler_angles(&q,order,intrinsic,&ang_x,&ang_y,&ang_z);
    exp_x = -1.4288992721907328e+00;
    exp_y = -3.3983690945412182e-01;
    exp_z = -2.3561944901923448e+00;
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    PDM_quaternion_to_euler_angles(&q,rev_order,extrinsic,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_to_rotation_matrix", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(1.,-2.,3.,-4,&q);
    PDM_quaternion_normalize(&q);
    double rot_mat[9];
    PDM_quaternion_to_rotation_matrix(&q,rot_mat);
    // printf("%s::%d\n [\n%23.16e,%23.16e,%23.16e,\n%23.16e,%23.16e,%23.16e,\n%23.16e,%23.16e,%23.16e]\n",
    //     __FILE__,__LINE__,
    //     rot_mat[3*0+0],
    //     rot_mat[3*0+1],
    //     rot_mat[3*0+2],
    //     rot_mat[3*1+0],
    //     rot_mat[3*1+1],
    //     rot_mat[3*1+2],
    //     rot_mat[3*2+0],
    //     rot_mat[3*2+1],
    //     rot_mat[3*2+2]);
    double exp_mat[9] = {
        -6.6666666666666663e-01,-1.3333333333333336e-01, 7.3333333333333317e-01,
        -6.6666666666666652e-01,-3.3333333333333331e-01,-6.6666666666666663e-01,
         3.3333333333333326e-01,-9.3333333333333324e-01, 1.3333333333333341e-01
    };
    CHECK_EQ_C_ARRAY_FLOAT(rot_mat,exp_mat,9,EPS);
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_to_rotation_matrix 2", 1) {
    PDM_quaternion q;
    PDM_quaternion q_out;
    PDM_quaternion_set(0.,0.,0.,0.,&q);
    PDM_quaternion_set(0.,0.,0.,0.,&q_out);
    PDM_bool_t intrinsic = PDM_TRUE;
    int order[3] = {2,1,0};
    double ang_x,ang_y,ang_z;
    double rot_mat[9];
    const double* r_mat = (const double*) rot_mat;

    for (int k = 0; k<3; k++) {
        ang_z = -10.+k*10.;
        for (int j = 0; j<3; j++) {
            ang_y = -10.+j*10.;
            for (int i = 0; i<3; i++) {
                ang_x = -10.+i*10.;
                PDM_quaternion_from_euler_angles(ang_x*DEG2RAD,ang_y*DEG2RAD,ang_z*DEG2RAD,order,intrinsic,&q);
                PDM_quaternion_to_rotation_matrix(&q,rot_mat);
                PDM_quaternion_from_rotation_matrix(r_mat,&q_out);
                CHECK(PDM_quaternion_equal_quaternion(&q,&q_out,EPS));
            }
        }
    }
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_to_rotation_matrix 3", 1) {
    PDM_quaternion q;
    PDM_quaternion_set(1.,-2.,3.,-4,&q);
    PDM_quaternion_normalize(&q);
    double rot_mat[9];
    PDM_quaternion_to_rotation_matrix(&q,rot_mat);
    // printf("%s::%d\n [\n%23.16e,%23.16e,%23.16e,\n%23.16e,%23.16e,%23.16e,\n%23.16e,%23.16e,%23.16e]\n",
    //     __FILE__,__LINE__,
    //     rot_mat[3*0+0],
    //     rot_mat[3*0+1],
    //     rot_mat[3*0+2],
    //     rot_mat[3*1+0],
    //     rot_mat[3*1+1],
    //     rot_mat[3*1+2],
    //     rot_mat[3*2+0],
    //     rot_mat[3*2+1],
    //     rot_mat[3*2+2]);
    double exp_mat[9] = {
        -6.6666666666666663e-01,-1.3333333333333336e-01, 7.3333333333333317e-01,
        -6.6666666666666652e-01,-3.3333333333333331e-01,-6.6666666666666663e-01,
         3.3333333333333326e-01,-9.3333333333333324e-01, 1.3333333333333341e-01
    };
    CHECK_EQ_C_ARRAY_FLOAT(rot_mat,exp_mat,9,EPS);
}
