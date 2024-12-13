#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_quaternion.h"

static const double EPS = 8*__DBL_EPSILON__;
static const double DEG2RAD = M_PI/180.;
static const double RAD2DEG = 180./M_PI;

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_set", 1) {
    double q[4];
    PDM_quaternion_set(0.,1.,2.,3.,q);
    CHECK(PDM_quaternion_equal(q,0.,1.,2.,3.,EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_equal", 1) {
    double q[4];
    PDM_quaternion_set(0.,1.,2.,3.,q);
    CHECK( PDM_quaternion_equal(q,0.,1.,2.,3.,EPS));
    CHECK(!PDM_quaternion_equal(q,1.,1.,2.,3.,EPS));
    CHECK(!PDM_quaternion_equal(q,0.,0.,2.,3.,EPS));
    CHECK(!PDM_quaternion_equal(q,0.,1.,0.,3.,EPS));
    CHECK(!PDM_quaternion_equal(q,0.,1.,2.,0.,EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_set_identity", 1) {
    double q[4];
    PDM_quaternion_set(0.,1.,2.,3.,q);
    PDM_quaternion_set_identity(q);
    CHECK(PDM_quaternion_equal(q,1.,0.,0.,0.,EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_conjugate", 1) {
    double q[4];
    double qt[4];
    PDM_quaternion_set(0.,1.,2.,3.,q);
    PDM_quaternion_set(3.,2.,1.,0.,qt);
    PDM_quaternion_conjugate(q,qt);
    CHECK(PDM_quaternion_equal(qt,q[0],-q[1],-q[2],-q[3],EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_norm", 1) {
    double q[4];
    PDM_quaternion_set(0.,3.,4.,0.,q);
    double norm = PDM_quaternion_norm(q);
    CHECK(PDM_ABS(norm-5.)<__DBL_EPSILON__);
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_normalize", 1) {
    double q[4];
    PDM_quaternion_set(0.,3.,4.,0.,q);
    PDM_quaternion_normalize(q);
    double norm = PDM_quaternion_norm(q);
    CHECK(PDM_ABS(norm-1.)<__DBL_EPSILON__);

    q[0] = 2.;
    q[1] = 0.;
    q[2] = 0.;
    q[3] = 0.;
    PDM_quaternion_normalize(q);
    CHECK(PDM_quaternion_equal(q,1.,0.,0.,0.,EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_rotate", 1) {
    // rotation of 120deg of axis 1,1,1 (i->k,k->j,j->i)
    double q[4];
    PDM_quaternion_set(-.5,.5,.5,.5,q);
    double vector[9] = {1.,0.,0.,
                        0.,1.,0.,
                        0.,0.,1.,
                        };
    double out_vect[9];
    double exp_vect[9] = {0.,0.,1.,
                          1.,0.,0.,
                          0.,1.,0.};
    PDM_quaternion_rotate(q,vector,3,out_vect);
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,9,EPS);

    double qt[4];
    PDM_quaternion_set(-.5,-.5,-.5,-.5,qt);
    double exp_vect_2[9] = {0.,1.,0.,
                            0.,0.,1.,
                            1.,0.,0.};
    PDM_quaternion_rotate(qt,vector,3,out_vect);
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect_2,9,EPS);
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_rotate 2", 1) {
    // rotation of 90 deg around y-axis
    double sin_45 = .5*sqrt(2.);
    double q[4];
    PDM_quaternion_set(sin_45,0.,sin_45,0.,q);
    double vector[3] = {1.,2.,3.};
    double out_vect[3];
    double exp_vect[3] = {3.,2.,-1.};
    PDM_quaternion_rotate(q,vector,1,out_vect);
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);

}

    
MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_compose", 1) {
    double q1[4];
    double q2[4];
    double q_out[4];
    PDM_quaternion_set(0.,0.,0.,0.,q1);
    PDM_quaternion_set(1.,2.,3.,4.,q2);
    PDM_quaternion_set(0.,0.,0.,0.,q_out);
    PDM_quaternion_set_identity(q1);
    PDM_quaternion_compose(q1,q2,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,q2,EPS));

    PDM_quaternion_compose(q2,q1,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,q2,EPS));

}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_compose 2", 1) {
    double q1[4];
    double qi[4];
    double qj[4];
    double qk[4];
    double q_out[4];
    PDM_quaternion_set(1.,0.,0.,0.,q1);
    PDM_quaternion_set(0.,1.,0.,0.,qi);
    PDM_quaternion_set(0.,0.,1.,0.,qj);
    PDM_quaternion_set(0.,0.,0.,1.,qk);
    PDM_quaternion_set(0.,0.,0.,0.,q_out);

    PDM_quaternion_compose(q1,q1,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,q1,EPS));
    PDM_quaternion_compose(q1,qi,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,qi,EPS));
    PDM_quaternion_compose(q1,qj,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,qj,EPS));
    PDM_quaternion_compose(q1,qk,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,qk,EPS));

    PDM_quaternion_compose(qi,q1,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,qi,EPS));
    PDM_quaternion_compose(qi,qi,q_out);
    CHECK(PDM_quaternion_equal(q_out,-q1[0],-q1[1],-q1[2],-q1[3],EPS));
    PDM_quaternion_compose(qi,qj,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,qk,EPS));
    PDM_quaternion_compose(qi,qk,q_out);
    CHECK(PDM_quaternion_equal(q_out,-qj[0],-qj[1],-qj[2],-qj[3],EPS));

    PDM_quaternion_compose(qj,q1,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,qj,EPS));
    PDM_quaternion_compose(qj,qi,q_out);
    CHECK(PDM_quaternion_equal(q_out,-qk[0],-qk[1],-qk[2],-qk[3],EPS));
    PDM_quaternion_compose(qj,qj,q_out);
    CHECK(PDM_quaternion_equal(q_out,-q1[0],-q1[1],-q1[2],-q1[3],EPS));
    PDM_quaternion_compose(qj,qk,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,qi,EPS));

    PDM_quaternion_compose(qk,q1,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,qk,EPS));
    PDM_quaternion_compose(qk,qi,q_out);
    CHECK(PDM_quaternion_equal_quaternion(q_out,qj,EPS));
    PDM_quaternion_compose(qk,qj,q_out);
    CHECK(PDM_quaternion_equal(q_out,-qi[0],-qi[1],-qi[2],-qi[3],EPS));
    PDM_quaternion_compose(qk,qk,q_out);
    CHECK(PDM_quaternion_equal(q_out,-q1[0],-q1[1],-q1[2],-q1[3],EPS));
    
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_two_vectors", 1) {
    double q[4];
    PDM_quaternion_set(0.,0.,0.,0.,q);
    double v1[3] = {1.,0.,0.};
    double v2[3] = {0.,1.,0.};
    PDM_quaternion_from_two_vectors(v1,v2,q);
    double sin_45 = .5*sqrt(2.);
    CHECK(PDM_quaternion_equal(q,sin_45,0.,0.,sin_45,EPS));

    double v3[3];
    PDM_quaternion_rotate(q,v1,1,v3);
    CHECK_EQ_C_ARRAY_FLOAT(v3,v2,3,EPS);
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_axis_angle", 1) {
    double q[4];
    PDM_quaternion_set(0.,0.,0.,0.,q);
    double axis[3] = {1.,1.,1.};
    double angle = 120.*M_PI/180.;
    PDM_quaternion_from_axis_angle(axis,angle,q);

    CHECK(PDM_quaternion_equal(q,.5,.5,.5,.5,EPS));

    double axis2[3] = {-1.,-1.,-1.};
    angle = -120.*M_PI/180.;
    PDM_quaternion_from_axis_angle(axis2,angle,q);

    CHECK(PDM_quaternion_equal(q,.5,.5,.5,.5,EPS));
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_axis_angle 2", 1) {
    double q[4];
    PDM_quaternion_set(0.,0.,0.,0.,q);
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
            PDM_quaternion_from_axis_angle(&axes[3*(5*i+j)],angle,q);
            CHECK(PDM_quaternion_equal(q,w[5*i+j],v0[5*i+j],v1[5*i+j],v2[5*i+j],EPS));
        }
    }
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_euler_angles", 1) {
    double q[4];
    PDM_quaternion_set(0.,0.,0.,0.,q);
    double sin_45 = .5*sqrt(2.);
    PDM_bool_t intrinsic = PDM_TRUE;
    int order[3] = {2,1,0};

    PDM_quaternion_from_euler_angles(.5*M_PI,0.,0.,order,intrinsic,q);
    CHECK(PDM_quaternion_equal(q,sin_45,sin_45,0.,0.,EPS));
    PDM_quaternion_from_euler_angles(0.,0.5*M_PI,0.,order,intrinsic,q);
    CHECK(PDM_quaternion_equal(q,sin_45,0.,sin_45,0.,EPS));
    PDM_quaternion_from_euler_angles(0.,0.,0.5*M_PI,order,intrinsic,q);
    CHECK(PDM_quaternion_equal(q,sin_45,0.,0.,sin_45,EPS));
    PDM_quaternion_from_euler_angles(0.5*M_PI,0.5*M_PI,0.,order,intrinsic,q);
    CHECK(PDM_quaternion_equal(q,.5,.5,.5,-.5,EPS));
    PDM_quaternion_from_euler_angles(0.5*M_PI,0.,0.5*M_PI,order,intrinsic,q);
    CHECK(PDM_quaternion_equal(q,.5,.5,.5,.5,EPS));
    PDM_quaternion_from_euler_angles(0.,0.5*M_PI,0.5*M_PI,order,intrinsic,q);
    CHECK(PDM_quaternion_equal(q,.5,-.5,.5,.5,EPS));

}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_euler_angles 2", 1) {
    double q[4];
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

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_from_rotation_matrix", 1) {
    double q[4];
    PDM_quaternion_set(0.,0.,0.,0.,q);
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

                PDM_quaternion_from_rotation_matrix(r_mat,q);
                PDM_quaternion_to_euler_angles(q,order,intrinsic,&out_ang_x,&out_ang_y,&out_ang_z);
                CHECK(PDM_ABS(ang_x-out_ang_x)<EPS);
                CHECK(PDM_ABS(ang_y-out_ang_y)<EPS);
                CHECK(PDM_ABS(ang_z-out_ang_z)<EPS);
            }
        }
    }
}

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_to_axis_angle", 1) {
    double q[4];
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
    double q[4];
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

MPI_TEST_CASE("[pdm_quaternion] - 1p - PDM_quaternion_to_rotation_matrix", 1) {
    double q[4];
    double q_out[4];
    PDM_quaternion_set(0.,0.,0.,0.,q);
    PDM_quaternion_set(0.,0.,0.,0.,q_out);
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
                PDM_quaternion_from_euler_angles(ang_x*DEG2RAD,ang_y*DEG2RAD,ang_z*DEG2RAD,order,intrinsic,q);
                PDM_quaternion_to_rotation_matrix(q,rot_mat);
                PDM_quaternion_from_rotation_matrix(r_mat,q_out);
                CHECK(PDM_quaternion_equal_quaternion(q,q_out,EPS));
            }
        }
    }
}

