#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_rotation.h"
/*----------------------------------------------------------------------------
 *  PDM_rotation_t INPUT/OUTPUT FUNCTIONS
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  PDM_rotation_t COMPOSITE FUNCTIONS
 *----------------------------------------------------------------------------*/

#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)


static const double EPS = 8*__DBL_EPSILON__;
static const double DEG2RAD = M_PI/180.;

static void multiply_matrices(double A[16],double B[16],double C[16]){
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            C[4*i+j] = 0;
            for (int k = 0; k < 4; k++) {
                C[4*i+j] += A[4*i+k] * B[4*k+j];
            }
        }
    }
}

static void mat_vec(double A[9],double x[3],double y[3]){
    for (int i = 0; i < 3; i++) {
        y[i] = 0;
        for (int j = 0; j < 3; j++) {
            y[i] += A[3*i+j]*x[j];
        }
    }
}

MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_multiply_n_by_n_matrices", 1) {
    double A[16] = {
         1.,  0.,  2.,  0.,
         0.,  1.,  0.,  0.,
        -2.,  0.,  1.,  0.,
         0.,  0.,  0.,  1.
    };
    double B[16] = {
        1., 0., 0., 0.,
        0., 1., 0., 0.,
        0., 0., 2., 0.,
        0., 0., 0., 1.
    };
    double D[16] = {
         1.,  2.,  3.,  4.,
         5.,  6.,  7.,  8.,
         9., 10., 11., 12.,
        13., 14., 15., 16.
    };
    double E[16] = {
        1.,-5., -9.,-13.,
        2., 6.,-10.,-14.,
        3., 7., 11.,-15.,
        4., 8., 12., 16.
    };
    double C[16];
    double exp_C[16];
    multiply_matrices(A,B,exp_C);
    PDM_rotation_multiply_n_by_n_matrices(A,B,4,C);
    CHECK_EQ_C_ARRAY_FLOAT(C,exp_C,12,EPS);

    multiply_matrices(D,E,exp_C);
    PDM_rotation_multiply_n_by_n_matrices(D,E,4,C);
    CHECK_EQ_C_ARRAY_FLOAT(C,exp_C,12,EPS);
}

MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_multiply_n_by_n_matrices_2", 1) {
    double A[16] = {
         1.,  2.,  3.,  4.,
         5.,  6.,  7.,  8.,
         9., 10., 11., 12.,
        13., 14., 15., 16.
    };
    double B[16] = {
        1.,-5., -9.,-13.,
        2., 6.,-10.,-14.,
        3., 7., 11.,-15.,
        4., 8., 12., 16.
    };
    double exp_C[16];
    double B_end[16];
    multiply_matrices(A,B,exp_C);
    PDM_rotation_multiply_n_by_n_matrices(A,B,4,B_end);
    CHECK_EQ_C_ARRAY_FLOAT(B_end,exp_C,12,EPS);
}

MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_multiply_n_by_n_matrices_3", 1) {
    double A[16] = {
         1.,  2.,  3.,  4.,
         5.,  6.,  7.,  8.,
         9., 10., 11., 12.,
        13., 14., 15., 16.
    };
    double B[16] = {
        1.,-5., -9.,-13.,
        2., 6.,-10.,-14.,
        3., 7., 11.,-15.,
        4., 8., 12., 16.
    };
    double exp_C[16];
    double A_end[16];
    multiply_matrices(A,B,exp_C);
    PDM_rotation_multiply_n_by_n_matrices(A,B,4,A_end);
    CHECK_EQ_C_ARRAY_FLOAT(A_end,exp_C,12,EPS);
}

MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_apply_n_by_n_matrix", 1) {

    double mat[9];
    double axis[3] = {1.,0.,1.};
    double angle = 12.*DEG2RAD;
    PDM_rotation_axis_angle_to_rotation_matrix(axis,angle,PDM_FALSE,mat);
    int n_samp = 7;
    double vector[21] = {
        1.,0.,0.,
        0.,1.,0.,
        0.,0.,1.,
        1.,1.,1.,
        0.,1.,1.,
        1.,0.,1.,
        1.,1.,0.
    };
    double vector_out[21];
    double exp_out[21];
    for (int i = 0; i < n_samp; i++){
        mat_vec(mat,&vector[3*i],&exp_out[3*i]);

    }
    PDM_rotation_apply_n_by_n_matrix(mat,vector,3,n_samp,vector_out);
    CHECK_EQ_C_ARRAY_FLOAT(vector_out,exp_out,12,EPS);
}

MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_axis_angle_to_euler_angles", 1) {
    double angle = 37*DEG2RAD;
    double axis[3] = {-3.,4.,-5};
    int order[3],rev_order[3];
    double ang_x,ang_y,ang_z;
    double exp_x,exp_y,exp_z;

    // direct, [0,1,2], intrinsic
    order[0] = 0;
    order[1] = 1;
    order[2] = 2;
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,order,PDM_TRUE,&ang_x,&ang_y,&ang_z);
    exp_x = -1.9195732121200226e-01;
    exp_y =  4.1244155008887562e-01;
    exp_z = -4.2437040657168013e-01;
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,rev_order,PDM_FALSE,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // direct, [0,2,1], intrinsic
    order[0] = 0;
    order[1] = 2;
    order[2] = 1;
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,order,PDM_TRUE,&ang_x,&ang_y,&ang_z);
    exp_x = -3.7112789994533268e-01;
    exp_y =  4.4762159939481483e-01;
    exp_z = -3.8679270256453857e-01;
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,rev_order,PDM_FALSE,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // direct, [1,0,2], intrinsic
    order[0] = 1;
    order[1] = 0;
    order[2] = 2;
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,order,PDM_TRUE,&ang_x,&ang_y,&ang_z);
    exp_x = -1.7568506091904279e-01;
    exp_y =  4.1929215641494677e-01;
    exp_z = -5.0211818050329748e-01;
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,rev_order,PDM_FALSE,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // direct, [1,2,0], intrinsic
    order[0] = 1;
    order[1] = 2;
    order[2] = 0;
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,order,PDM_TRUE,&ang_x,&ang_y,&ang_z);
    exp_x = -1.9981002701339715e-01;
    exp_y =  3.2361964206307259e-01;
    exp_z = -4.9368599869889129e-01;
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,rev_order,PDM_FALSE,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // direct, [2,0,1], intrinsic
    order[0] = 2;
    order[1] = 0;
    order[2] = 1;
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,order,PDM_TRUE,&ang_x,&ang_y,&ang_z);
    exp_x = -3.4253328559234930e-01;
    exp_y =  3.0186325235851291e-01;
    exp_z = -4.1204660950517663e-01;
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,rev_order,PDM_FALSE,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // direct, [2,1,0], intrinsic
    order[0] = 2;
    order[1] = 1;
    order[2] = 0;
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,order,PDM_TRUE,&ang_x,&ang_y,&ang_z);
    exp_x = -3.5743456350361774e-01;
    exp_y =  2.8382394280353651e-01;
    exp_z = -5.1625197447159565e-01;
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
    for (int i=0;i<3;i++){
        rev_order[i] = order[2-i];
    }
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,rev_order,PDM_FALSE,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // reverse
    PDM_rotation_axis_angle_to_euler_angles(axis,-angle,PDM_FALSE,order,PDM_TRUE,&exp_x,&exp_y,&exp_z);
    PDM_rotation_axis_angle_to_euler_angles(axis, angle,PDM_TRUE, order,PDM_TRUE,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_axis_angle_to_rotation_matrix", 1) {
    double angle = 37*DEG2RAD;
    double axis[3] = {-3.,4.,-5};
    double rot_mat[9];
    double exp_mat[9];
    // testing only reverse
    PDM_rotation_axis_angle_to_rotation_matrix(axis,-angle,PDM_FALSE,rot_mat);
    PDM_rotation_axis_angle_to_rotation_matrix(axis, angle,PDM_TRUE, exp_mat);
    CHECK_EQ_C_ARRAY_FLOAT(rot_mat,exp_mat,9,EPS);
}

MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_axis_angle_to_homogeneous_matrix", 1) {
    double angle = 37*DEG2RAD;
    double axis[3] = {-3.,4.,-5};
    double homo_mat[16];
    double expe_mat[16];
    // testing only reverse
    PDM_rotation_axis_angle_to_homogeneous_matrix(axis,-angle,PDM_FALSE,homo_mat);
    PDM_rotation_axis_angle_to_homogeneous_matrix(axis, angle,PDM_TRUE, expe_mat);
    CHECK_EQ_C_ARRAY_FLOAT(homo_mat,expe_mat,16,EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix", 1) {
    double angle = 37*DEG2RAD;
    double axis[3] = {-3.,4.,-5};
    double homo_mat[16];
    double expe_mat[16];
    double rot_center[3] = {1,-2,3};
    double orig[3] = {0,0,0};
    PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix(axis,angle,orig,PDM_FALSE,expe_mat);
    PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix(axis,angle,rot_center,PDM_FALSE,homo_mat);
    double rot_mat[9];
    double exp_mat[9];
    PDM_rotation_homogeneous_matrix_to_rotation_matrix(homo_mat,PDM_FALSE,rot_mat);
    PDM_rotation_homogeneous_matrix_to_rotation_matrix(expe_mat,PDM_FALSE,exp_mat);
    CHECK_EQ_C_ARRAY_FLOAT(rot_mat,exp_mat,9,EPS);
    double RT[3];
    PDM_rotation_apply_n_by_n_matrix(rot_mat,rot_center,3,1,RT);
    double expec_T[3] = {
        rot_center[0]-RT[0],
        rot_center[1]-RT[1],
        rot_center[2]-RT[2],
    };
    CHECK(PDM_ABS(homo_mat[3] -expec_T[0]) < EPS);
    CHECK(PDM_ABS(homo_mat[7] -expec_T[1]) < EPS);
    CHECK(PDM_ABS(homo_mat[11]-expec_T[2]) < EPS);

    // reverse
    PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix(axis,angle,rot_center,PDM_FALSE,expe_mat);
    PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix(axis,-angle,rot_center,PDM_TRUE,homo_mat);
    CHECK_EQ_C_ARRAY_FLOAT(homo_mat,expe_mat,16,EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_euler_angles_to_axis_angle", 1) {
    double ang_x = 15*DEG2RAD;
    double ang_y = -25*DEG2RAD;
    double ang_z = 135*DEG2RAD;
    int order[3];
    PDM_bool_t intrinsic = PDM_TRUE;
    double axis[3];
    double angle;
    double expec_ax[3];
    double expec_angle;

    // direct, [0,1,2], intrinsic
    order[0] = 0;
    order[1] = 1;
    order[2] = 2;
    PDM_rotation_euler_angles_to_axis_angle(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,axis,&angle);
    expec_ax[0] = -1.6283521498290401e-01;
    expec_ax[1] = -2.1769635115701913e-01;
    expec_ax[2] =  9.6233725452898844e-01;
    expec_angle =  2.3261541784711866e+00;
    CHECK_EQ_C_ARRAY_FLOAT(axis,expec_ax,3,EPS);
    CHECK(PDM_ABS(angle-expec_angle)<EPS);

    // direct, [1,0,2], intrinsic
    order[0] = 1;
    order[1] = 0;
    order[2] = 2;
    PDM_rotation_euler_angles_to_axis_angle(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,axis,&angle);
    expec_ax[0] = -1.5922306689701238e-01;
    expec_ax[1] = -2.1286722707461483e-01;
    expec_ax[2] =  9.6402051773054798e-01;
    expec_angle =  2.4385735552120056e+00;
    CHECK_EQ_C_ARRAY_FLOAT(axis,expec_ax,3,EPS);
    CHECK(PDM_ABS(angle-expec_angle)<EPS);

    // direct, [2,1,0], intrinsic
    order[0] = 2;
    order[1] = 1;
    order[2] = 0;
    PDM_rotation_euler_angles_to_axis_angle(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,axis,&angle);
    expec_ax[0] =  2.6310757008103358e-01;
    expec_ax[1] =  3.7932149432640026e-02;
    expec_ax[2] =  9.6402051773054809e-01;
    expec_angle =  2.4385735552120051e+00;
    CHECK_EQ_C_ARRAY_FLOAT(axis,expec_ax,3,EPS);
    CHECK(PDM_ABS(angle-expec_angle)<EPS);

    // reverse, [2,1,0], intrinsic
    PDM_rotation_euler_angles_to_axis_angle(ang_x,ang_y,ang_z,order,intrinsic,PDM_TRUE,axis,&angle);
    expec_ax[0] = -2.6310757008103358e-01;
    expec_ax[1] = -3.7932149432640026e-02;
    expec_ax[2] = -9.6402051773054809e-01;
    expec_angle =  2.4385735552120051e+00;
    CHECK_EQ_C_ARRAY_FLOAT(axis,expec_ax,3,EPS);
    CHECK(PDM_ABS(angle-expec_angle)<EPS);
}

MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_euler_angles_to_euler_angles", 1) {
    double ang_x = 15*DEG2RAD;
    double ang_y = -25*DEG2RAD;
    double ang_z = 135*DEG2RAD;
    int order[3];
    int out_order[3];
    double out_ang_x,out_ang_y,out_ang_z;
    double exp_ang_x,exp_ang_y,exp_ang_z;
    PDM_bool_t intrinsic = PDM_TRUE;

    order[0] = 0;
    order[1] = 1;
    order[2] = 2;
    // direct, [0,1,2], intrinsic => [0,1,2], intrinsic
    out_order[0] = 0;
    out_order[1] = 1;
    out_order[2] = 2;
    PDM_rotation_euler_angles_to_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,out_order,intrinsic,&out_ang_x,&out_ang_y,&out_ang_z);
    exp_ang_x = 2.6179938779914957e-01;
    exp_ang_y =-4.3633231299858233e-01;
    exp_ang_z = 2.3561944901923448e+00;
    CHECK(PDM_ABS(ang_x-exp_ang_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_ang_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_ang_z)<EPS);

    // direct, [0,1,2], intrinsic => [1,0,2], intrinsic
    out_order[0] = 1;
    out_order[1] = 0;
    out_order[2] = 2;
    PDM_rotation_euler_angles_to_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,out_order,intrinsic,&out_ang_x,&out_ang_y,&out_ang_z);
    exp_ang_x = 2.3677591065618109e-01;
    exp_ang_y =-4.4975847000749192e-01;
    exp_ang_z = 2.2434346174088571e+00;
    CHECK(PDM_ABS(out_ang_x-exp_ang_x)<EPS);
    CHECK(PDM_ABS(out_ang_y-exp_ang_y)<EPS);
    CHECK(PDM_ABS(out_ang_z-exp_ang_z)<EPS);

    // direct, [0,1,2], intrinsic => [2,1,0], intrinsic
    out_order[0] = 2;
    out_order[1] = 1;
    out_order[2] = 0;
    PDM_rotation_euler_angles_to_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,out_order,intrinsic,&out_ang_x,&out_ang_y,&out_ang_z);
    exp_ang_x =-4.9419199674160330e-01;
    exp_ang_y = 1.0583842460540160e-01;
    exp_ang_z = 2.2711165511661227e+00;
    CHECK(PDM_ABS(out_ang_x-exp_ang_x)<EPS);
    CHECK(PDM_ABS(out_ang_y-exp_ang_y)<EPS);
    CHECK(PDM_ABS(out_ang_z-exp_ang_z)<EPS);

    order[0] = 2;
    order[1] = 1;
    order[2] = 0;
    // direct, [2,1,0], intrinsic => [0,2,1], intrinsic
    out_order[0] = 0;
    out_order[1] = 2;
    out_order[2] = 1;
    PDM_rotation_euler_angles_to_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,out_order,intrinsic,&out_ang_x,&out_ang_y,&out_ang_z);
    exp_ang_x = 2.8423566328092438e+00;
    exp_ang_y = 2.5071159211902962e+00;
    exp_ang_z = 6.5060533465895043e-01;
    CHECK(PDM_ABS(out_ang_x-exp_ang_x)<EPS);
    CHECK(PDM_ABS(out_ang_y-exp_ang_y)<EPS);
    CHECK(PDM_ABS(out_ang_z-exp_ang_z)<EPS);

}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_euler_angles_to_rotation_matrix", 1) {
    double ang_x = 15*DEG2RAD;
    double ang_y = -25*DEG2RAD;
    double ang_z = 135*DEG2RAD;
    int order[3];
    order[0] = 2;
    order[1] = 1;
    order[2] = 0;
    PDM_bool_t intrinsic = PDM_TRUE;

    double rot_mat[9];
    double exp_mat[9]= {
        -6.4085638205578854e-01, -6.0566819194206101e-01,  4.7166634271272756e-01,
         6.4085638205578865e-01, -7.6035721184237781e-01, -1.0564093892828903e-01,
         4.2261826174069950e-01,  2.3456971600980447e-01,  8.7542609806559313e-01
    };
    PDM_rotation_euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,rot_mat);
    CHECK_EQ_C_ARRAY_FLOAT(rot_mat,exp_mat,9,EPS);

    PDM_rotation_euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z,order,intrinsic,PDM_TRUE,rot_mat);
    double exp_mat_rev[9] = {
        -6.4085638205578854e-01,  6.4085638205578865e-01,  4.2261826174069950e-01,
        -6.0566819194206101e-01, -7.6035721184237781e-01,  2.3456971600980447e-01,
         4.7166634271272756e-01, -1.0564093892828903e-01,  8.7542609806559313e-01
    };
    CHECK_EQ_C_ARRAY_FLOAT(rot_mat,exp_mat_rev,9,EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_euler_angles_to_homogeneous_matrix", 1) {
    double ang_x = 15*DEG2RAD;
    double ang_y = -25*DEG2RAD;
    double ang_z = 135*DEG2RAD;
    int order[3];
    order[0] = 2;
    order[1] = 1;
    order[2] = 0;
    PDM_bool_t intrinsic = PDM_TRUE;

    double rot_mat[16];
    double exp_mat[16]= {
        -6.4085638205578854e-01, -6.0566819194206101e-01,  4.7166634271272756e-01, 0.,
         6.4085638205578865e-01, -7.6035721184237781e-01, -1.0564093892828903e-01, 0.,
         4.2261826174069950e-01,  2.3456971600980447e-01,  8.7542609806559313e-01, 0.,
         0.,0.,0.,1.
    };
    PDM_rotation_euler_angles_to_homogeneous_matrix(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,rot_mat);
    CHECK_EQ_C_ARRAY_FLOAT(rot_mat,exp_mat,16,EPS);

    PDM_rotation_euler_angles_to_homogeneous_matrix(ang_x,ang_y,ang_z,order,intrinsic,PDM_TRUE,rot_mat);
    double exp_mat_rev[16] = {
        -6.4085638205578854e-01,  6.4085638205578865e-01,  4.2261826174069950e-01, 0.,
        -6.0566819194206101e-01, -7.6035721184237781e-01,  2.3456971600980447e-01, 0.,
         4.7166634271272756e-01, -1.0564093892828903e-01,  8.7542609806559313e-01, 0.,
         0.,0.,0.,1.
    };
    CHECK_EQ_C_ARRAY_FLOAT(rot_mat,exp_mat_rev,16,EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_euler_angles_and_rotation_center_to_homogeneous_matrix", 1) {
    double ang_x = 15*DEG2RAD;
    double ang_y = -25*DEG2RAD;
    double ang_z = 135*DEG2RAD;
    int order[3];
    order[0] = 2;
    order[1] = 1;
    order[2] = 0;
    PDM_bool_t intrinsic = PDM_TRUE;
    double homo_mat[16];
    double rot_mat[9];
    double exp_mat[9];
    double rot_center[3] = {1,-2,3};
    PDM_rotation_euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,exp_mat);
    PDM_rotation_euler_angles_and_rotation_center_to_homogeneous_matrix(ang_x,ang_y,ang_z,order,intrinsic,rot_center,PDM_FALSE,homo_mat);
    PDM_rotation_homogeneous_matrix_to_rotation_matrix(homo_mat,PDM_FALSE,rot_mat);
    CHECK_EQ_C_ARRAY_FLOAT(rot_mat,exp_mat,9,EPS);
    double RT[3];
    PDM_rotation_apply_n_by_n_matrix(rot_mat,rot_center,3,1,RT);
    double expec_T[3] = {
        rot_center[0]-RT[0],
        rot_center[1]-RT[1],
        rot_center[2]-RT[2],
    };
    CHECK(PDM_ABS(homo_mat[3] -expec_T[0]) < EPS);
    CHECK(PDM_ABS(homo_mat[7] -expec_T[1]) < EPS);
    CHECK(PDM_ABS(homo_mat[11]-expec_T[2]) < EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_periodic_t_info_to_homogeneous_matrix", 1) {
    double rotation_center[3] = { 1.,-2., 3.};
    double rotation_angle[3] = {.4,-.5,.6};
    double translation[3] = {-4, 5,-6};
    double hmat[16];
    double rmat[9];
    double rot_mat[9];
    int order[3] = {2,1,0};
    double out_center[3],out_angle[3],out_trsl[3];

    PDM_rotation_periodic_t_info_to_homogeneous_matrix(rotation_center,rotation_angle,translation,PDM_FALSE,hmat);
    PDM_rotation_euler_angles_to_rotation_matrix(rotation_angle[0],rotation_angle[1],rotation_angle[2],order,PDM_TRUE,PDM_FALSE,rot_mat);
    PDM_rotation_homogeneous_matrix_to_rotation_matrix(hmat,PDM_FALSE,rmat);
    PDM_rotation_homogeneous_matrix_to_periodic_t_info(hmat,PDM_FALSE,PDM_FALSE,out_center,out_angle,out_trsl);
    CHECK_EQ_C_ARRAY_FLOAT(rmat,rot_mat,9,EPS);

    CHECK(PDM_ABS(hmat[3]  - out_trsl[0])<EPS);
    CHECK(PDM_ABS(hmat[7]  - out_trsl[1])<EPS);
    CHECK(PDM_ABS(hmat[11] - out_trsl[2])<EPS);

    PDM_rotation_periodic_t_info_to_homogeneous_matrix(rotation_center,rotation_angle,translation,PDM_TRUE,hmat);
    PDM_rotation_euler_angles_to_rotation_matrix(rotation_angle[0],rotation_angle[1],rotation_angle[2],order,PDM_TRUE,PDM_TRUE,rot_mat);
    PDM_rotation_homogeneous_matrix_to_rotation_matrix(hmat,PDM_FALSE,rmat);
    PDM_rotation_homogeneous_matrix_to_periodic_t_info(hmat,PDM_FALSE,PDM_FALSE,out_center,out_angle,out_trsl);
    CHECK_EQ_C_ARRAY_FLOAT(rmat,rot_mat,9,EPS);

    CHECK(PDM_ABS(hmat[3] - out_trsl[0])<EPS);
    CHECK(PDM_ABS(hmat[7] - out_trsl[1])<EPS);
    CHECK(PDM_ABS(hmat[11] - out_trsl[2])<EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_rotation_matrix_to_axis_angle", 1) {
    double axis[3] = {1.,-2.,3.};
    double angle = 135*DEG2RAD;
    double rot_mat[9];
    double out_axis[3];
    double exp_axis[3];
    double out_ang;
    PDM_rotation_axis_angle_to_rotation_matrix(axis,angle,PDM_FALSE,rot_mat);

    PDM_rotation_rotation_matrix_to_axis_angle(rot_mat,PDM_FALSE,out_axis,&out_ang);
    exp_axis[0] =  2.6726124191242445e-01;
    exp_axis[1] = -5.3452248382484879e-01;
    exp_axis[2] =  8.0178372573727308e-01;

    CHECK_EQ_C_ARRAY_FLOAT(out_axis,exp_axis,3,EPS);
    CHECK(PDM_ABS(out_ang-angle)<EPS);

    PDM_rotation_rotation_matrix_to_axis_angle(rot_mat,PDM_TRUE,out_axis,&out_ang);
    exp_axis[0] = -2.6726124191242445e-01;
    exp_axis[1] =  5.3452248382484879e-01;
    exp_axis[2] = -8.0178372573727308e-01;

    CHECK_EQ_C_ARRAY_FLOAT(out_axis,exp_axis,3,EPS);
    CHECK(PDM_ABS(out_ang-angle)<EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_rotation_matrix_to_euler_angles", 1) {
    double axis[3] = {1.,-2.,3.};
    double angle = 135*DEG2RAD;
    double rot_mat[9];
    double ang_x,ang_y,ang_z;
    double exp_x,exp_y,exp_z;
    PDM_rotation_axis_angle_to_rotation_matrix(axis,angle,PDM_FALSE,rot_mat);


    PDM_bool_t intrinsic = PDM_TRUE;
    int order[3];

    // [0,2,1], intrinsic
    order[0] = 0;
    order[1] = 2;
    order[2] = 1;
    PDM_rotation_rotation_matrix_to_euler_angles(rot_mat,PDM_FALSE,order,intrinsic,
                                                &ang_x,&ang_y,&ang_z);
    exp_x = -1.9549639796314939e+00;
    exp_y = -3.1208224212373898e+00;
    exp_z =  9.4555023416975281e-01;
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // [0,2,1], intrinsic reverse
    PDM_rotation_rotation_matrix_to_euler_angles(rot_mat,PDM_TRUE,order,intrinsic,
                                                &ang_x,&ang_y,&ang_z);
    exp_x = -1.8047159482796502e+00;
    exp_y =  2.2374115338684435e+00;
    exp_z = -3.2897621227233076e-01;
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // [2,1,0], intrinsic
    order[0] = 2;
    order[1] = 1;
    order[2] = 0;
    PDM_rotation_rotation_matrix_to_euler_angles(rot_mat,PDM_FALSE,order,intrinsic,
                                                &ang_x,&ang_y,&ang_z);
    exp_x = -9.4723239412904581e-01;
    exp_y = -8.3869742630145583e-01;
    exp_z =  2.6371364489318241e+00;
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);

    // [2,1,0], intrinsic reverse
    PDM_rotation_rotation_matrix_to_euler_angles(rot_mat,PDM_TRUE,order,intrinsic,
                                                &ang_x,&ang_y,&ang_z);
    exp_x = -1.1697869430919412e+00;
    exp_y =  1.2156176430153876e-02;
    exp_z = -2.1959400502139310e+00;
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_rotation_matrix_to_homogeneous_matrix", 1) {
    double axis[3] = {1.,-2.,3.};
    double angle = 135*DEG2RAD;
    double rot_mat[9];
    double homo_mat[16];
    PDM_rotation_axis_angle_to_rotation_matrix(axis,angle,PDM_FALSE,rot_mat);
    PDM_rotation_rotation_matrix_to_homogeneous_matrix(rot_mat,PDM_FALSE,homo_mat);
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            CHECK(PDM_ABS(rot_mat[3*i+j]-homo_mat[4*i+j])<EPS);
        }
    }
    CHECK(PDM_ABS(homo_mat[4*0+3])<EPS);
    CHECK(PDM_ABS(homo_mat[4*1+3])<EPS);
    CHECK(PDM_ABS(homo_mat[4*2+3])<EPS);
    CHECK(PDM_ABS(homo_mat[4*3+3]-1)<EPS);
    CHECK(PDM_ABS(homo_mat[4*3+0])<EPS);
    CHECK(PDM_ABS(homo_mat[4*3+1])<EPS);
    CHECK(PDM_ABS(homo_mat[4*3+2])<EPS);

    PDM_rotation_rotation_matrix_to_homogeneous_matrix(rot_mat,PDM_TRUE,homo_mat);
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            CHECK(PDM_ABS(rot_mat[3*i+j]-homo_mat[4*j+i])<EPS);
        }
    }
    CHECK(PDM_ABS(homo_mat[4*0+3])<EPS);
    CHECK(PDM_ABS(homo_mat[4*1+3])<EPS);
    CHECK(PDM_ABS(homo_mat[4*2+3])<EPS);
    CHECK(PDM_ABS(homo_mat[4*3+3]-1)<EPS);
    CHECK(PDM_ABS(homo_mat[4*3+0])<EPS);
    CHECK(PDM_ABS(homo_mat[4*3+1])<EPS);
    CHECK(PDM_ABS(homo_mat[4*3+2])<EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix", 1) {
    double axis[3] = {1.,-2.,3.};
    double angle = 135*DEG2RAD;
    double rot_mat[9];
    double rot_center[3] = {-4.,5.,-6.};
    double homo_mat[16];
    PDM_rotation_axis_angle_to_rotation_matrix(axis,angle,PDM_FALSE,rot_mat);

    PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix(rot_mat,rot_center,PDM_FALSE,homo_mat);
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            CHECK(PDM_ABS(rot_mat[3*i+j]-homo_mat[4*i+j])<EPS);
        }
    }
    double RT[3];
    PDM_rotation_apply_n_by_n_matrix(rot_mat,rot_center,3,1,RT);
    double expec_T[3] = {
        rot_center[0]-RT[0],
        rot_center[1]-RT[1],
        rot_center[2]-RT[2],
    };
    CHECK(PDM_ABS(homo_mat[3] -expec_T[0]) < EPS);
    CHECK(PDM_ABS(homo_mat[7] -expec_T[1]) < EPS);
    CHECK(PDM_ABS(homo_mat[11]-expec_T[2]) < EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_homogeneous_matrix_to_axis_angle", 1) {
    double axis[3] = {1.,-2.,3.};
    double angle = 135*DEG2RAD;
    double homo_mat[16];

    double out_axis[3];
    double exp_axis[3];
    double out_ang;

    PDM_rotation_axis_angle_to_homogeneous_matrix(axis,angle,PDM_FALSE,homo_mat);
    PDM_rotation_homogeneous_matrix_to_axis_angle(homo_mat,PDM_FALSE,out_axis,&out_ang);
    exp_axis[0] =  2.6726124191242445e-01;
    exp_axis[1] = -5.3452248382484879e-01;
    exp_axis[2] =  8.0178372573727308e-01;

    CHECK_EQ_C_ARRAY_FLOAT(out_axis,exp_axis,3,EPS);
    CHECK(PDM_ABS(out_ang-angle)<EPS);

    PDM_rotation_homogeneous_matrix_to_axis_angle(homo_mat,PDM_TRUE,out_axis,&out_ang);
    exp_axis[0] = -2.6726124191242445e-01;
    exp_axis[1] =  5.3452248382484879e-01;
    exp_axis[2] = -8.0178372573727308e-01;

    CHECK_EQ_C_ARRAY_FLOAT(out_axis,exp_axis,3,EPS);
    CHECK(PDM_ABS(out_ang-angle)<EPS);

}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_homogeneous_matrix_to_euler_angles", 1) {
    double axis[3] = {1.,-2.,3.};
    double angle = 135*DEG2RAD;
    double rot_mat[9];
    double homo_mat[16];
    PDM_bool_t intrinsic = PDM_TRUE;
    int order[3];
    double ang_x,ang_y,ang_z;
    double exp_x,exp_y,exp_z;
    PDM_rotation_axis_angle_to_rotation_matrix(axis,angle,PDM_FALSE,rot_mat);
    PDM_rotation_axis_angle_to_homogeneous_matrix(axis,angle,PDM_FALSE,homo_mat);

    // [0,1,2]
    order[0] = 0;
    order[1] = 1;
    order[2] = 2;
    PDM_rotation_rotation_matrix_to_euler_angles(rot_mat,PDM_FALSE,order,intrinsic,&exp_x,&exp_y,&exp_z);
    PDM_rotation_homogeneous_matrix_to_euler_angles(homo_mat,PDM_FALSE,order,intrinsic,&ang_x,&ang_y,&ang_z);
    CHECK(PDM_ABS(ang_x-exp_x)<EPS);
    CHECK(PDM_ABS(ang_y-exp_y)<EPS);
    CHECK(PDM_ABS(ang_z-exp_z)<EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_homogeneous_matrix_to_euler_angles_and_translation", 1) {
    double axis[3] = {1.,-2.,3.};
    double angle = 135*DEG2RAD;
    double homo_mat[16];
    double translation[3];
    double exp_trsl[3];
    int order[3];
    PDM_bool_t intrinsic = PDM_TRUE;
    double ang_x,ang_y,ang_z;
    PDM_rotation_axis_angle_to_homogeneous_matrix(axis,angle,PDM_FALSE,homo_mat);
    homo_mat[3]  = -4.;
    homo_mat[7]  =  5.;
    homo_mat[11] = -6.;

    // [2,1,0]
    order[0] = 2;
    order[1] = 1;
    order[2] = 0;
    PDM_rotation_homogeneous_matrix_to_euler_angles_and_translation(homo_mat,PDM_FALSE,order,intrinsic,&ang_x,&ang_y,&ang_z,translation);
    exp_trsl[0] = -4.0000000000000000e+00;
    exp_trsl[1] =  5.0000000000000000e+00;
    exp_trsl[2] = -6.0000000000000000e+00;
    CHECK_EQ_C_ARRAY_FLOAT(translation,exp_trsl,3,EPS);
    PDM_rotation_homogeneous_matrix_to_euler_angles_and_translation(homo_mat,PDM_TRUE,order,intrinsic,&ang_x,&ang_y,&ang_z,translation);
    exp_trsl[0] =  5.0658452273779386e-01;
    exp_trsl[1] = -5.4022762270905922e+00;
    exp_trsl[2] =  6.8962876743603436e+00;
    CHECK_EQ_C_ARRAY_FLOAT(translation,exp_trsl,3,10*EPS); // sensitive to compilation options
}

MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_homogeneous_matrix_to_periodic_t_info no rot_center", 1) {
    double axis[3] = {1.,-2.,3.};
    double angle = 135*DEG2RAD;
    double homo_mat[16];
    double ang_x,ang_y,ang_z;
    int order[3] = {2,1,0};
    PDM_bool_t intrinsic = PDM_TRUE;
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,order,intrinsic,&ang_x,&ang_y,&ang_z);
    PDM_rotation_axis_angle_to_homogeneous_matrix(axis,angle,PDM_FALSE,homo_mat);
    homo_mat[3]  = -4.;
    homo_mat[7]  =  5.;
    homo_mat[11] = -6.;
    PDM_bool_t compute_rotation_center = PDM_FALSE;

    double rotation_center[3];
    double rotation_angle[3];
    double translation[3];

    PDM_rotation_homogeneous_matrix_to_periodic_t_info(homo_mat,PDM_FALSE,compute_rotation_center,
                                                       rotation_center,rotation_angle,translation);
    CHECK(PDM_ABS(rotation_center[0])<EPS);
    CHECK(PDM_ABS(rotation_center[1])<EPS);
    CHECK(PDM_ABS(rotation_center[2])<EPS);
    CHECK(PDM_ABS(rotation_angle[0]-ang_x)<EPS);
    CHECK(PDM_ABS(rotation_angle[1]-ang_y)<EPS);
    CHECK(PDM_ABS(rotation_angle[2]-ang_z)<EPS);
    CHECK(PDM_ABS(translation[0]-homo_mat[3])<EPS);
    CHECK(PDM_ABS(translation[1]-homo_mat[7])<EPS);
    CHECK(PDM_ABS(translation[2]-homo_mat[11])<EPS);

    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_TRUE,order,intrinsic,&ang_x,&ang_y,&ang_z);
    PDM_rotation_homogeneous_matrix_to_periodic_t_info(homo_mat,PDM_TRUE,compute_rotation_center,
                                                       rotation_center,rotation_angle,translation);
    CHECK(PDM_ABS(rotation_center[0])<EPS);
    CHECK(PDM_ABS(rotation_center[1])<EPS);
    CHECK(PDM_ABS(rotation_center[2])<EPS);
    CHECK(PDM_ABS(rotation_angle[0]-ang_x)<EPS);
    CHECK(PDM_ABS(rotation_angle[1]-ang_y)<EPS);
    CHECK(PDM_ABS(rotation_angle[2]-ang_z)<EPS);
    double exp_trsl[3];
    exp_trsl[0] =  5.0658452273779386e-01;
    exp_trsl[1] = -5.4022762270905922e+00;
    exp_trsl[2] =  6.8962876743603436e+00;
    CHECK_EQ_C_ARRAY_FLOAT(translation,exp_trsl,3,10*EPS); // sensitive to compilation options
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_homogeneous_matrix_to_periodic_t_info rot_center", 1) {
    double axis[3] = {1.,-2.,3.};
    double angle = 135*DEG2RAD;
    double homo_mat[16];
    double rot_mat[9];
    double ang_x,ang_y,ang_z;
    int order[3] = {2,1,0};
    PDM_bool_t intrinsic = PDM_TRUE;
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_FALSE,order,intrinsic,&ang_x,&ang_y,&ang_z);
    PDM_rotation_axis_angle_to_homogeneous_matrix(axis,angle,PDM_FALSE,homo_mat);
    PDM_rotation_axis_angle_to_rotation_matrix(axis,angle,PDM_FALSE,rot_mat);
    homo_mat[3]  = -4.;
    homo_mat[7]  =  5.;
    homo_mat[11] = -6.;
    PDM_bool_t compute_rotation_center = PDM_TRUE;

    double rotation_center[3];
    double rotation_angle[3];
    double translation[3];
    double cross[3];
    double check[3];

    PDM_rotation_homogeneous_matrix_to_periodic_t_info(homo_mat,PDM_FALSE,compute_rotation_center,
                                                       rotation_center,rotation_angle,translation);
    PDM_rotation_apply_n_by_n_matrix(rot_mat,rotation_center,3,1,check);

    // checking that -RC+R+T == hmat[:,-1]
    check[0] *= -1;
    check[0] += rotation_center[0]+translation[0];
    check[1] *= -1;
    check[1] += rotation_center[1]+translation[1];
    check[2] *= -1;
    check[2] += rotation_center[2]+translation[2];
    CHECK(PDM_ABS(check[0]-homo_mat[3])<EPS);
    CHECK(PDM_ABS(check[1]-homo_mat[7])<EPS);
    CHECK(PDM_ABS(check[2]-homo_mat[11])<EPS);
    // checking translation is // to axis
    PDM_CROSS_PRODUCT(cross,axis,translation); // should be 0
    CHECK(PDM_ABS(cross[0])<2*EPS);
    CHECK(PDM_ABS(cross[1])<2*EPS);
    CHECK(PDM_ABS(cross[2])<2*EPS);
    CHECK(PDM_ABS(rotation_angle[0]-ang_x)<EPS);
    CHECK(PDM_ABS(rotation_angle[1]-ang_y)<EPS);
    CHECK(PDM_ABS(rotation_angle[2]-ang_z)<EPS);
    double exp_trsl[3];
    exp_trsl[0] = -2.2857142857142865e+00;
    exp_trsl[1] =  4.5714285714285712e+00;
    exp_trsl[2] = -6.8571428571428577e+00;
    CHECK_EQ_C_ARRAY_FLOAT(translation,exp_trsl,3,EPS);

    // reverse
    PDM_rotation_axis_angle_to_euler_angles(axis,angle,PDM_TRUE,order,intrinsic,&ang_x,&ang_y,&ang_z);
    PDM_rotation_homogeneous_matrix_to_periodic_t_info(homo_mat,PDM_TRUE,compute_rotation_center,
                                                    rotation_center,rotation_angle,translation);
    PDM_rotation_apply_n_by_n_matrix(rot_mat,rotation_center,3,1,check);
    // checking that -RC+R+T == hmat[:,-1]
    check[0] *= -1;
    check[0] += rotation_center[0]-translation[0];
    check[1] *= -1;
    check[1] += rotation_center[1]-translation[1];
    check[2] *= -1;
    check[2] += rotation_center[2]-translation[2];
    CHECK(PDM_ABS(check[0]-homo_mat[3])<EPS);
    CHECK(PDM_ABS(check[1]-homo_mat[7])<EPS);
    CHECK(PDM_ABS(check[2]-homo_mat[11])<EPS);
    // checking translation is // to axis
    PDM_CROSS_PRODUCT(cross,axis,translation); // should be 0
    CHECK(PDM_ABS(cross[0])<2*EPS);
    CHECK(PDM_ABS(cross[1])<2*EPS);
    CHECK(PDM_ABS(cross[2])<2*EPS);
    CHECK(PDM_ABS(rotation_angle[0]-ang_x)<EPS);
    CHECK(PDM_ABS(rotation_angle[1]-ang_y)<EPS);
    CHECK(PDM_ABS(rotation_angle[2]-ang_z)<EPS);
    exp_trsl[0] =  2.2857142857142865e+00;
    exp_trsl[1] = -4.5714285714285712e+00;
    exp_trsl[2] =  6.8571428571428577e+00;
    CHECK_EQ_C_ARRAY_FLOAT(translation,exp_trsl,3,EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_homogeneous_matrix_to_rotation_matrix", 1) {
    double axis[3] = {1.,-2.,3.};
    double angle = 135*DEG2RAD;
    double homo_mat[16];
    double rot_mat[9];
    PDM_rotation_axis_angle_to_homogeneous_matrix(axis,angle,PDM_FALSE,homo_mat);
    PDM_rotation_homogeneous_matrix_to_rotation_matrix(homo_mat,PDM_FALSE,rot_mat);
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            CHECK(PDM_ABS(rot_mat[3*i+j]-homo_mat[4*i+j])<EPS);
        }
    }
    PDM_rotation_homogeneous_matrix_to_rotation_matrix(homo_mat,PDM_TRUE,rot_mat);
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            CHECK(PDM_ABS(rot_mat[3*i+j]-homo_mat[4*j+i])<EPS);
        }
    }
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_two_vectors_to_axis_angle", 1) {
    double vect_1[3] = {-1., 2.,-3.};
    double vect_2[3] = { 4.,-5., 6.};
    double axis[3];
    double ang;
    double out_vect[3];
    double exp_vect[3];
    double rot_mat[9];

    PDM_rotation_two_vectors_to_axis_angle(vect_1,vect_2,PDM_FALSE,axis,&ang);
    PDM_rotation_axis_angle_to_rotation_matrix(axis,ang,PDM_FALSE,rot_mat);
    PDM_rotation_apply_n_by_n_matrix(rot_mat,vect_1,3,1,out_vect);
    exp_vect[0] =  1.7056057308448840e+00;
    exp_vect[1] = -2.1320071635561049e+00;
    exp_vect[2] =  2.5584085962673235e+00;
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);

    PDM_rotation_two_vectors_to_axis_angle(vect_1,vect_2,PDM_TRUE,axis,&ang);
    PDM_rotation_axis_angle_to_rotation_matrix(axis,ang,PDM_FALSE,rot_mat);
    PDM_rotation_apply_n_by_n_matrix(rot_mat,vect_2,3,1,out_vect);
    exp_vect[0] = -2.3452078799117118e+00;
    exp_vect[1] =  4.6904157598234306e+00;
    exp_vect[2] = -7.0356236397351442e+00;
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,10*EPS); // sensitive to compilation options

    // permuting vect_1 & vect_2
    double axis_2[3];
    double ang_2;
    PDM_rotation_two_vectors_to_axis_angle(vect_2,vect_1,PDM_FALSE,axis_2,&ang_2);
    CHECK_EQ_C_ARRAY_FLOAT(axis_2,axis,3,EPS);
    CHECK(PDM_ABS(ang-ang_2)<EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_two_vectors_to_euler_angles", 1) {
    double vect_1[3] = {-1., 2.,-3.};
    double vect_2[3] = { 4.,-5., 6.};
    double ang_x,ang_y,ang_z;
    int order[3] = {2,1,0};
    PDM_bool_t intrinsic = PDM_TRUE;
    double out_vect[3];
    double exp_vect[3];
    double rot_mat[9];

    PDM_rotation_two_vectors_to_euler_angles(vect_1,vect_2,PDM_FALSE,order,intrinsic,&ang_x,&ang_y,&ang_z);
    PDM_rotation_euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,rot_mat);
    PDM_rotation_apply_n_by_n_matrix(rot_mat,vect_1,3,1,out_vect);
    exp_vect[0] =  1.7056057308448838e+00;
    exp_vect[1] = -2.1320071635561044e+00;
    exp_vect[2] =  2.5584085962673262e+00;
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);

    // reverse
    PDM_rotation_two_vectors_to_euler_angles(vect_1,vect_2,PDM_TRUE,order,intrinsic,&ang_x,&ang_y,&ang_z);
    PDM_rotation_euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,rot_mat);
    PDM_rotation_apply_n_by_n_matrix(rot_mat,vect_2,3,1,out_vect);
    exp_vect[0] = -2.3452078799117126e+00;
    exp_vect[1] =  4.6904157598234297e+00;
    exp_vect[2] = -7.0356236397351450e+00;
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);

}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_two_vectors_to_rotation_matrix", 1) {
    double vect_1[3] = {-1., 2.,-3.};
    double vect_2[3] = { 4.,-5., 6.};
    double rot_mat[9];
    double out_vect[3];
    double exp_vect[3];
    PDM_rotation_two_vectors_to_rotation_matrix(vect_1,vect_2,PDM_FALSE,rot_mat);
    PDM_rotation_apply_n_by_n_matrix(rot_mat,vect_1,3,1,out_vect);
    exp_vect[0] =  1.7056057308448838e+00;
    exp_vect[1] = -2.1320071635561044e+00;
    exp_vect[2] =  2.5584085962673262e+00;
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);

    // reverse
    PDM_rotation_two_vectors_to_rotation_matrix(vect_1,vect_2,PDM_TRUE,rot_mat);
    PDM_rotation_apply_n_by_n_matrix(rot_mat,vect_2,3,1,out_vect);
    exp_vect[0] = -2.3452078799117153e+00;
    exp_vect[1] =  4.6904157598234315e+00;
    exp_vect[2] = -7.0356236397351477e+00;
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);

}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_two_vectors_to_homogeneous_matrix", 1) {
    double vect_1[3] = {-1., 2.,-3.};
    double vect_2[3] = { 4.,-5., 6.};
    double homo_mat[16];
    double out_vect[3];
    double exp_vect[3];
    PDM_rotation_two_vectors_to_homogeneous_matrix(vect_1,vect_2,PDM_FALSE,homo_mat);
    PDM_rotation_apply_homogeneous_matrix(homo_mat,vect_1,1,out_vect);
    exp_vect[0] =  1.7056057308448838e+00;
    exp_vect[1] = -2.1320071635561044e+00;
    exp_vect[2] =  2.5584085962673262e+00;
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);

    // reverse
    PDM_rotation_two_vectors_to_homogeneous_matrix(vect_1,vect_2,PDM_TRUE,homo_mat);
    PDM_rotation_apply_homogeneous_matrix(homo_mat,vect_2,1,out_vect);
    exp_vect[0] = -2.3452078799117153e+00;
    exp_vect[1] =  4.6904157598234315e+00;
    exp_vect[2] = -7.0356236397351477e+00;
    CHECK_EQ_C_ARRAY_FLOAT(out_vect,exp_vect,3,EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_two_vectors_and_rotation_center_to_homogeneous_matrix", 1) {
    double vect_1[3] = {-1., 2.,-3.};
    double vect_2[3] = { 4.,-5., 6.};
    double homo_mat[16];
    double rot_mat[9];
    double exp_mat[9];
    double rot_center[3] = {1,-2,3};
    PDM_rotation_two_vectors_to_rotation_matrix(vect_1,vect_2,PDM_FALSE,exp_mat);
    PDM_rotation_two_vectors_and_rotation_center_to_homogeneous_matrix(vect_1,vect_2,rot_center,PDM_FALSE,homo_mat);
    PDM_rotation_homogeneous_matrix_to_rotation_matrix(homo_mat,PDM_FALSE,rot_mat);
    CHECK_EQ_C_ARRAY_FLOAT(rot_mat,exp_mat,9,EPS);
    double RT[3];
    PDM_rotation_apply_n_by_n_matrix(rot_mat,rot_center,3,1,RT);
    double expec_T[3] = {
        rot_center[0]-RT[0],
        rot_center[1]-RT[1],
        rot_center[2]-RT[2],
    };
    CHECK(PDM_ABS(homo_mat[3] -expec_T[0]) < EPS);
    CHECK(PDM_ABS(homo_mat[7] -expec_T[1]) < EPS);
    CHECK(PDM_ABS(homo_mat[11]-expec_T[2]) < EPS);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_axes_and_origin_to_homogeneous_matrix", 1) {
    double angle = 37*DEG2RAD;
    double axis[3] = {-3.,4.,-5};
    double rot_mat[9];
    double e[9] = {
        1.,0.,0.,
        0.,1.,0.,
        0.,0.,1.,
    };
    PDM_rotation_axis_angle_to_rotation_matrix(axis,angle,PDM_FALSE,rot_mat);
    double r[9];
    double in[3];
    double out[3];
    PDM_rotation_apply_n_by_n_matrix(rot_mat,e,3,3,r);
    double axis_1[3] = {r[3*0+0],r[3*0+1],r[3*0+2]};
    double axis_2[3] = {r[3*1+0],r[3*1+1],r[3*1+2]};
    double axis_3[3] = {r[3*2+0],r[3*2+1],r[3*2+2]};
    double orig[3] = {-4., 5,-6.};

    double hmat[16];
    PDM_rotation_axes_and_origin_to_homogeneous_matrix(axis_1,axis_2,axis_3,orig,PDM_FALSE,hmat);

    // testing conversion of B base vectors in A system to B system
    in[0] = axis_1[0]+orig[0];
    in[1] = axis_1[1]+orig[1];
    in[2] = axis_1[2]+orig[2];
    PDM_rotation_apply_homogeneous_matrix(hmat,in,1,out);
    CHECK_EQ_C_ARRAY_FLOAT(&e[0],out,3,2*EPS);

    in[0] = axis_2[0]+orig[0];
    in[1] = axis_2[1]+orig[1];
    in[2] = axis_2[2]+orig[2];
    PDM_rotation_apply_homogeneous_matrix(hmat,in,1,out);
    CHECK_EQ_C_ARRAY_FLOAT(&e[3],out,3,2*EPS);

    in[0] = axis_3[0]+orig[0];
    in[1] = axis_3[1]+orig[1];
    in[2] = axis_3[2]+orig[2];
    PDM_rotation_apply_homogeneous_matrix(hmat,in,1,out);
    CHECK_EQ_C_ARRAY_FLOAT(&e[6],out,3,2*EPS);

    double hmat_rev[16];
    PDM_rotation_axes_and_origin_to_homogeneous_matrix(axis_1,axis_2,axis_3,orig,PDM_TRUE,hmat_rev);
    // testing conversion of B base vectors in B system to A system
    in[0] = 1.;
    in[1] = 0.;
    in[2] = 0.;
    PDM_rotation_apply_homogeneous_matrix(hmat_rev,in,1,out);
    out[0] -= orig[0];
    out[1] -= orig[1];
    out[2] -= orig[2];
    CHECK_EQ_C_ARRAY_FLOAT(out,axis_1,3,2*EPS);

    in[0] = 0.;
    in[1] = 1.;
    in[2] = 0.;
    PDM_rotation_apply_homogeneous_matrix(hmat_rev,in,1,out);
    out[0] -= orig[0];
    out[1] -= orig[1];
    out[2] -= orig[2];
    CHECK_EQ_C_ARRAY_FLOAT(out,axis_2,3,2*EPS);

    in[0] = 0.;
    in[1] = 0.;
    in[2] = 1.;
    PDM_rotation_apply_homogeneous_matrix(hmat_rev,in,1,out);
    out[0] -= orig[0];
    out[1] -= orig[1];
    out[2] -= orig[2];
    CHECK_EQ_C_ARRAY_FLOAT(out,axis_3,3,2*EPS);

}

MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_apply_euler_angles_and_rotation_center", 1) {
    PDM_bool_t intrinsic = PDM_TRUE;
    PDM_bool_t reverse = PDM_FALSE;
    int order[3] = {2,1,0};
    double ang_x = 5.*DEG2RAD;
    double ang_y = 10.*DEG2RAD;
    double ang_z = -15.*DEG2RAD;
    double rotation_center[3] = {1.,2.,3.};
    int n_samp = 4;
    double vector[12] = {
        1.,0.,0.,
        0.,1.,0.,
        0.,0.,1.,
        1.,1.,1.
    };
    double expec_vector_out[12];
    double tmp_vec[3];
    double tmp_mat[9];
    // doing it by hand:
    PDM_rotation_euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z,order,intrinsic,PDM_FALSE,tmp_mat);
    // print_matrix(tmp_mat,3,3);
    for (int i = 0; i < n_samp; i++){
        // applying the translation
        for (int j = 0; j < 3; j++){
            tmp_vec[j] = vector[3*i+j] - rotation_center[j];
        }
        // applying the rotation
        mat_vec(tmp_mat,tmp_vec,&expec_vector_out[3*i]);
        // applying the translation back
        for (int j = 0; j < 3; j++){
            expec_vector_out[3*i+j] += rotation_center[j];
        }
    }

    double vector_out[12];
    PDM_rotation_apply_euler_angles_and_rotation_center(ang_x,ang_y,ang_z,order,
        intrinsic,rotation_center,reverse,vector,n_samp,vector_out);

    CHECK_EQ_C_ARRAY_FLOAT(vector_out,expec_vector_out,12,EPS);
}

MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_apply_axis_angle_and_rotation_center", 1) {
    double axis[3] = {1.,2.,3.};
    double angle = 12.*DEG2RAD;
    double rotation_center[3] = {-1.,2.,-3.};
    PDM_bool_t reverse = PDM_FALSE;
    int n_samp = 4;
    double vector[12] = {
        1.,0.,0.,
        0.,1.,0.,
        0.,0.,1.,
        1.,1.,1.
    };
    double expec_vector_out[12];

    double tmp_vec[3];
    double tmp_mat[9];
    double vector_out[12];
    // doing it by hand:
    PDM_rotation_axis_angle_to_rotation_matrix(axis,angle,PDM_FALSE,tmp_mat);
    // print_matrix(tmp_mat,3,3);
    for (int i = 0; i < n_samp; i++){
        // applying the translation
        for (int j = 0; j < 3; j++){
            tmp_vec[j] = vector[3*i+j] - rotation_center[j];
        }
        // applying the rotation
        mat_vec(tmp_mat,tmp_vec,&expec_vector_out[3*i]);
        // applying the translation back
        for (int j = 0; j < 3; j++){
            expec_vector_out[3*i+j] += rotation_center[j];
        }
    }
    PDM_rotation_apply_axis_angle_and_rotation_center(axis,angle,rotation_center,
        reverse,vector,n_samp,vector_out);
    CHECK_EQ_C_ARRAY_FLOAT(vector_out,expec_vector_out,12,EPS);
}

MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_apply_rotation_matrix_and_rotation_center", 1) {
    double axis[3] = {1.,2.,3.};
    double angle = 12.*DEG2RAD;
    double rotation_center[3] = {-1.,2.,-3.};
    PDM_bool_t reverse = PDM_FALSE;
    int n_samp = 4;
    double vector[12] = {
        1.,0.,0.,
        0.,1.,0.,
        0.,0.,1.,
        1.,1.,1.
    };
    double expec_vector_out[12];

    double tmp_vec[3];
    double tmp_mat[9];
    double vector_out[12];
    // doing it by hand:
    PDM_rotation_axis_angle_to_rotation_matrix(axis,angle,PDM_FALSE,tmp_mat);
    // print_matrix(tmp_mat,3,3);
    for (int i = 0; i < n_samp; i++){
        // applying the translation
        for (int j = 0; j < 3; j++){
            tmp_vec[j] = vector[3*i+j] - rotation_center[j];
        }
        // applying the rotation
        mat_vec(tmp_mat,tmp_vec,&expec_vector_out[3*i]);
        // applying the translation back
        for (int j = 0; j < 3; j++){
            expec_vector_out[3*i+j] += rotation_center[j];
        }
    }
    PDM_rotation_apply_rotation_matrix_and_rotation_center(tmp_mat,rotation_center,
        reverse,vector,n_samp,vector_out);
    CHECK_EQ_C_ARRAY_FLOAT(vector_out,expec_vector_out,12,EPS);
}

#endif
