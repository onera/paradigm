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
static const double RAD2DEG = 180./M_PI;

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

#endif
