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

// static void print_matrix(const double* mat,const int n_row,const int n_col) {
//     PDM_printf("#####\n");
//     for (int i=0; i<n_row;i++){
//         PDM_printf("|");
//         for (int j=0; j<n_col;j++){
//             PDM_printf("%5.3f\t",mat[n_col*i+j]);
//         }
//         PDM_printf("|\n");
//     }
//     PDM_printf("#####\n");
// }

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
    // print_matrix(C,4,4);
    // print_matrix(exp_C,4,4);
    CHECK_EQ_C_ARRAY_FLOAT(C,exp_C,12,EPS);

    multiply_matrices(D,E,exp_C);
    PDM_rotation_multiply_n_by_n_matrices(D,E,4,C);
    // print_matrix(C,4,4);
    // print_matrix(exp_C,4,4);
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
    // print_matrix(mat,3,3);
    // print_matrix(vector,7,3);
    // print_matrix(vector_out,7,3);
    // print_matrix(exp_out,7,3);
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
    // printf("%s::%d intrinsic order[%d %d %d] => [%23.16e %23.16e %23.16e]\n",__FILE__,__LINE__,
    //     order[0],order[1],order[2],ang_x,ang_y,ang_z);
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
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_euler_angles_to_euler_angles", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_euler_angles_to_rotation_matrix", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_euler_angles_to_homogeneous_matrix", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_euler_angles_and_rotation_center_to_homogeneous_matrix", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_periodic_t_info_to_homogeneous_matrix", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_rotation_matrix_to_axis_angle", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_rotation_matrix_to_euler_angles", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_rotation_matrix_to_homogeneous_matrix", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_homogeneous_matrix_to_axis_angle", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_homogeneous_matrix_to_euler_angles", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_homogeneous_matrix_to_euler_angles_and_translation", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_homogeneous_matrix_to_periodic_t_info", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_homogeneous_matrix_to_rotation_matrix", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_two_vectors_to_axis_angle", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_two_vectors_to_euler_angles", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_two_vectors_to_rotation_matrix", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_two_vectors_to_homogeneous_matrix", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_two_vectors_and_rotation_center_to_homogeneous_matrix", 1) {
    CHECK(0==1);
}
MPI_TEST_CASE("[pdm_rotation] - 1p - PDM_rotation_axes_and_origin_to_homogeneous_matrix", 1) {
    CHECK(0==1);
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
    CHECK(0==1);
}

#endif
