/*============================================================================
 * Mesure des temps CPU et elapsed
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <sysexits.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
/* after pdm.h: pdm.h includes pdm_config.h which defines PDM_HAVE_MKL and PDM_HAVE_LAPACK */
#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
#include <cblas.h>
#endif
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_quaternion.h"
#include "pdm_rotation.h"


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/


static
void 
set_identity_to_homogeneous_matrix
(
  double* homogeneous_matrix
)
{
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 4; j++){
      if (i==j){
        homogeneous_matrix[4*i+j] = 1.;
      }else{
        homogeneous_matrix[4*i+j] = 0.;
      }
    }
  }
}

// static 
// void 
// print_matrix(const double* mat,const int n_row,const int n_col) {
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

static
void
set_translation_to_homogeneous_matrix
(
  const double translation_vector[3],
  PDM_bool_t reverse,
        double homogeneous_matrix[16]
)
{
  set_identity_to_homogeneous_matrix(homogeneous_matrix);
  if (reverse){
    homogeneous_matrix[3]  = -translation_vector[0];
    homogeneous_matrix[7]  = -translation_vector[1];
    homogeneous_matrix[11] = -translation_vector[2];
  }else{
    homogeneous_matrix[3]  = translation_vector[0];
    homogeneous_matrix[7]  = translation_vector[1];
    homogeneous_matrix[11] = translation_vector[2];
  }
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------
 *  PDM_rotation_t INPUT/OUTPUT FUNCTIONS
 *----------------------------------------------------------------------------*/

void 
PDM_rotation_axis_angle_to_euler_angles
(
  const double axis[3],
  const double angle,
  const int order[3],
  const PDM_bool_t intrinsic,
        double* ang_x,
        double* ang_y,
        double* ang_z
)
{
  double qt[4];
  PDM_quaternion_from_axis_angle(axis,angle,qt);
  PDM_quaternion_to_euler_angles(qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void 
PDM_rotation_axis_angle_to_rotation_matrix
(
  const double axis[3],
  const double angle,
        double *rotation_matrix
)
{
  double qt[4];
  PDM_quaternion_from_axis_angle(axis,angle,qt);
  PDM_quaternion_to_rotation_matrix(qt,rotation_matrix);
}

void 
PDM_rotation_axis_angle_to_homogeneous_matrix
(
  const double axis[3],
  const double angle,
        double *homogeneous_matrix
)
{
  double qt[4];
  PDM_quaternion_from_axis_angle(axis,angle,qt);
  PDM_quaternion_to_homogeneous_matrix(qt,homogeneous_matrix);
}

void
PDM_rotation_euler_angles_to_axis_angle
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  PDM_bool_t intrinsic,
        double axis[3],
        double* angle
)
{
  double qt[4];
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,qt);
  PDM_quaternion_to_axis_angle(qt,axis,angle);
}

void
PDM_rotation_euler_angles_to_euler_angles
(
  const double input_ang_x,
  const double input_ang_y,
  const double input_ang_z,
  const int input_order[3],
  const PDM_bool_t input_intrinsic,
  const int output_order[3],
  const PDM_bool_t output_intrinsic,
        double* output_ang_x,
        double* output_ang_y,
        double* output_ang_z
)
{
  double qt[4];
  PDM_quaternion_from_euler_angles(input_ang_x,input_ang_y,input_ang_z,input_order,input_intrinsic,qt);
  PDM_quaternion_to_euler_angles(qt,output_order,output_intrinsic,output_ang_x,output_ang_y,output_ang_z);
}

void
PDM_rotation_euler_angles_to_rotation_matrix
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  PDM_bool_t intrinsic,
        double* rotation_matrix
)
{
  double qt[4];
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,qt);
  PDM_quaternion_to_rotation_matrix(qt,rotation_matrix);
}

void
PDM_rotation_euler_angles_to_homogeneous_matrix
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  PDM_bool_t intrinsic,
        double* homogeneous_matrix
)
{
  double qt[4];
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,qt);
  PDM_quaternion_to_homogeneous_matrix(qt,homogeneous_matrix);
}


void 
PDM_rotation_rotation_matrix_to_axis_angle
(
  const double* rotation_matrix,
        double axis[3],
        double* angle
)
{
  double qt[4];
  PDM_quaternion_from_rotation_matrix(rotation_matrix,qt);
  PDM_quaternion_to_axis_angle(qt,axis,angle);

}

void 
PDM_rotation_rotation_matrix_to_euler_angles
(
  const double* rotation_matrix,
  const int order[3],
  const PDM_bool_t intrinsic,
        double* ang_x,
        double* ang_y,
        double* ang_z
)
{
  double qt[4];
  PDM_quaternion_from_rotation_matrix(rotation_matrix,qt);
  PDM_quaternion_to_euler_angles(qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void 
PDM_rotation_rotation_matrix_to_homogeneous_matrix
(
  const double* rotation_matrix,
        double* homogeneous_matrix
)
{
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      homogeneous_matrix[4*i+j] = rotation_matrix[3*i+j];
    }
  }
  homogeneous_matrix[4*0+3] = 0.;
  homogeneous_matrix[4*1+3] = 0.;
  homogeneous_matrix[4*2+3] = 0.;
  homogeneous_matrix[4*3+0] = 0.;
  homogeneous_matrix[4*3+1] = 0.;
  homogeneous_matrix[4*3+2] = 0.;
  homogeneous_matrix[4*3+3] = 1.;
}

void 
PDM_rotation_homogeneous_matrix_to_axis_angle
(
  const double* homogeneous_matrix,
        double axis[3],
        double* angle
)
{
  double qt[4];
  PDM_quaternion_from_homogeneous_matrix(homogeneous_matrix,qt);
  PDM_quaternion_to_axis_angle(qt,axis,angle);
}

void 
PDM_rotation_homogeneous_matrix_to_euler_angles
(
  const double* homogeneous_matrix,
  const int order[3],
  const PDM_bool_t intrinsic,
        double* ang_x,
        double* ang_y,
        double* ang_z
)
{
  double qt[4];
  PDM_quaternion_from_homogeneous_matrix(homogeneous_matrix,qt);
  PDM_quaternion_to_euler_angles(qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void 
PDM_rotation_homogeneous_matrix_to_rotation_matrix
(
  const double* homogeneous_matrix,
        double* rotation_matrix
)
{
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      rotation_matrix[3*i+j] = homogeneous_matrix[4*i+j];
    }
  }
}

void 
PDM_rotation_two_vectors_to_axis_angle
(
  const double vector_1[3],
  const double vector_2[3],
        double axis[3],
        double* angle
)
{
  double qt[4];
  PDM_quaternion_from_two_vectors(vector_1,vector_2,qt);
  PDM_quaternion_to_axis_angle(qt,axis,angle);
}

void 
PDM_rotation_two_vectors_to_euler_angles
(
  const double vector_1[3],
  const double vector_2[3],
  const int order[3],
  const PDM_bool_t intrinsic,
        double* ang_x,
        double* ang_y,
        double* ang_z
)
{
  double qt[4];
  PDM_quaternion_from_two_vectors(vector_1,vector_2,qt);
  PDM_quaternion_to_euler_angles(qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void 
PDM_rotation_two_vectors_to_rotation_matrix
(
  const double vector_1[3],
  const double vector_2[3],
        double *rotation_matrix
)
{
  double qt[4];
  PDM_quaternion_from_two_vectors(vector_1,vector_2,qt);
  PDM_quaternion_to_rotation_matrix(qt,rotation_matrix);
}

void 
PDM_rotation_two_vectors_to_homogeneous_matrix
(
  const double vector_1[3],
  const double vector_2[3],
        double *homogeneous_matrix
)
{
  double qt[4];
  PDM_quaternion_from_two_vectors(vector_1,vector_2,qt);
  PDM_quaternion_to_homogeneous_matrix(qt,homogeneous_matrix);
}

/*----------------------------------------------------------------------------
 *  PDM_rotation_t HOMOGENOUS MATRICES
 *----------------------------------------------------------------------------*/
void
PDM_rotation_multiply_n_by_n_matrices
(
  const double* A,
  const double* B,
  const int     n,
        double* C
)
{
#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
  cblas_dgemm(CblasRowMajor,
              CblasNoTrans,
              CblasNoTrans,
              n,
              n,
              n,
              1.,
              A,
              n,
              B,
              n,
              0.,
              C,
              n);
#else
  PDM_UNUSED(A);
  PDM_UNUSED(B);
  PDM_UNUSED(n);
  PDM_UNUSED(C);
  printf("Error : CBLAS is mandatory (shipped with LAPACK or MKL), recompile with it.\n");
  exit(EX_CONFIG);
#endif
}

void 
PDM_rotation_apply_n_by_n_matrix
(
  const double *A,
  const double *x,
  const int     n,
  const int     n_samp,
        double *y
)
{
#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
  cblas_dgemm(CblasRowMajor,
              CblasNoTrans,
              CblasTrans,
              n_samp,
              n,n,
              1.,
              x,n,
              A,
              n,
              0.,
              y,
              n);
#else
  PDM_UNUSED(A);
  PDM_UNUSED(x);
  PDM_UNUSED(n);
  PDM_UNUSED(n_samp);
  PDM_UNUSED(y);
  printf("Error : CBLAS is mandatory (shipped with LAPACK or MKL), recompile with it.\n");
  exit(EX_CONFIG);
#endif
}

void 
PDM_rotation_apply_homogeneous_matrix
(
  const double homogeneous_matrix[16],
  const double* vector,
  const int n_samp,
        double* vector_out
)
{
  // applying the rotation
  double rotation_matrix[9] = {
    homogeneous_matrix[0], homogeneous_matrix[1], homogeneous_matrix[2],
    homogeneous_matrix[4], homogeneous_matrix[5], homogeneous_matrix[6],
    homogeneous_matrix[8], homogeneous_matrix[9], homogeneous_matrix[10],
  };
  PDM_rotation_apply_n_by_n_matrix(rotation_matrix,vector,3,n_samp,vector_out);

  // applying the translation
  for (int i = 0; i < n_samp; i++) {
    vector_out[3*i+0] += homogeneous_matrix[3];
    vector_out[3*i+1] += homogeneous_matrix[7];
    vector_out[3*i+2] += homogeneous_matrix[11];
  }
}

void
PDM_rotation_compose_homogeneous_matrices
(
  const double** homogeneous_matrices,
  const int      n_matrices,
        double   output_matrix[16]
)
{
  set_identity_to_homogeneous_matrix(output_matrix);
  for (int i = 0; i < n_matrices; i++) {
    double output_matrix_tmp[16];
    PDM_rotation_multiply_n_by_n_matrices(homogeneous_matrices[i],
                                          output_matrix,
                                          4,
                                          output_matrix_tmp);
    for(int j = 0; j < 16; ++j) {
      output_matrix[j] = output_matrix_tmp[j];
    }
  }
}

/*----------------------------------------------------------------------------
 *  PDM_rotation_t COMPOSITE FUNCTIONS
 *----------------------------------------------------------------------------*/

void 
PDM_rotation_apply_euler_angles_and_rotation_center
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  const PDM_bool_t intrinsic,
  const double rotation_center[3],
  const PDM_bool_t reverse,
  const double* vector,
  const int n_samp,
        double* vector_out
)
{
  // building the homogeneous matrix corresponding to whole transformation:
  // Trans+.Rot.Trans-
  double homogeneous_matrix    [16];
  double homogeneous_matrix_tmp[16];
  double tmp_matrix[16];
  // translation of -rotation_center
  set_translation_to_homogeneous_matrix(rotation_center,PDM_TRUE,
    homogeneous_matrix);
  
  // rotation
  PDM_rotation_euler_angles_to_homogeneous_matrix(ang_x,ang_y,ang_z,order,
    intrinsic,tmp_matrix);
  if (reverse){
    // transposing the rotation matrix
    double tmp;
    tmp = tmp_matrix[1];
    tmp_matrix[1] = tmp_matrix[4];
    tmp_matrix[4] = tmp;
    tmp = tmp_matrix[2];
    tmp_matrix[2] = tmp_matrix[8];
    tmp_matrix[8] = tmp;
    tmp = tmp_matrix[6];
    tmp_matrix[6] = tmp_matrix[9];
    tmp_matrix[9] = tmp;
  }
  PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
                                        homogeneous_matrix,
                                        4,
                                        homogeneous_matrix_tmp);
  for(int i = 0; i < 16; ++i) {
    homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  }
  
  // translation of rotation_center
  set_translation_to_homogeneous_matrix(rotation_center,
                                        PDM_FALSE,
                                        tmp_matrix);
  PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
                                        homogeneous_matrix,
                                        4,
                                        homogeneous_matrix_tmp);
  for(int i = 0; i < 16; ++i) {
    homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  }

  // applying the homogeneous matrix to the coordinate vector
  PDM_rotation_apply_homogeneous_matrix(homogeneous_matrix,
                                        vector,
                                        n_samp,
                                        vector_out);

}

void
PDM_rotation_apply_axis_angle_and_rotation_center
(
  const double axis[3],
  const double angle,
  const double rotation_center[3],
  const PDM_bool_t reverse,
  const double* vector,
  const int n_samp,
        double* vector_out
)
{
  // building the homogeneous matrix corresponding to whole transformation:
  // Trans+.Rot.Trans-
  double homogeneous_matrix[16];
  double homogeneous_matrix_tmp[16];
  double tmp_matrix[16];
  // translation of -rotation_center
  set_translation_to_homogeneous_matrix(rotation_center,PDM_TRUE,
    homogeneous_matrix);
  
  // rotation
  PDM_rotation_axis_angle_to_homogeneous_matrix(axis,angle,tmp_matrix);
  if (reverse){
    // transposing the rotation matrix
    double tmp;
    tmp = tmp_matrix[1];
    tmp_matrix[1] = tmp_matrix[4];
    tmp_matrix[4] = tmp;
    tmp = tmp_matrix[2];
    tmp_matrix[2] = tmp_matrix[8];
    tmp_matrix[8] = tmp;
    tmp = tmp_matrix[6];
    tmp_matrix[6] = tmp_matrix[9];
    tmp_matrix[9] = tmp;
  }

  // BLAS_DGEMM does not support inplace
  PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
                                        homogeneous_matrix,
                                        4,
                                        homogeneous_matrix_tmp);

  for(int i = 0; i < 16; ++i) {
    homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  }
  
  // translation of rotation_center
  set_translation_to_homogeneous_matrix(rotation_center,
                                        PDM_FALSE,
                                        tmp_matrix);

  PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
                                        homogeneous_matrix,
                                        4,
                                        homogeneous_matrix_tmp);

  for(int i = 0; i < 16; ++i) {
    homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  }

  // applying the homogeneous matrix to the coordinate vector
  PDM_rotation_apply_homogeneous_matrix(homogeneous_matrix,vector,n_samp,
    vector_out);
}

void
PDM_rotation_apply_rotation_matrix_and_rotation_center
(
  const double rotation_matrix[9],
  const double rotation_center[3],
  const PDM_bool_t reverse,
  const double* vector,
  const int n_samp,
        double* vector_out
)
{
  // building the homogeneous matrix corresponding to whole transformation:
  // Trans+.Rot.Trans-
  double homogeneous_matrix[16];
  double homogeneous_matrix_tmp[16];
  double tmp_matrix[16];
  // translation of -rotation_center
  set_translation_to_homogeneous_matrix(rotation_center,PDM_TRUE,
    homogeneous_matrix);
  
  // rotation
  PDM_rotation_rotation_matrix_to_homogeneous_matrix(rotation_matrix,
    tmp_matrix);
  if (reverse){
    // transposing the rotation matrix
    double tmp;
    tmp = tmp_matrix[1];
    tmp_matrix[1] = tmp_matrix[4];
    tmp_matrix[4] = tmp;
    tmp = tmp_matrix[2];
    tmp_matrix[2] = tmp_matrix[8];
    tmp_matrix[8] = tmp;
    tmp = tmp_matrix[6];
    tmp_matrix[6] = tmp_matrix[9];
    tmp_matrix[9] = tmp;
  }
  PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
                                        homogeneous_matrix,
                                        4,
                                        homogeneous_matrix_tmp);
  
  for(int i = 0; i < 16; ++i) {
    homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  }

  // translation of rotation_center
  set_translation_to_homogeneous_matrix(rotation_center,
                                        PDM_FALSE,
                                        tmp_matrix);

  PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
                                        homogeneous_matrix,
                                        4,
                                        homogeneous_matrix_tmp);

  for(int i = 0; i < 16; ++i) {
    homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  }

  // applying the homogeneous matrix to the coordinate vector
  PDM_rotation_apply_homogeneous_matrix(homogeneous_matrix,vector,n_samp,
    vector_out);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
