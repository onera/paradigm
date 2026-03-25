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
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_quaternion.h"
#include "pdm_quaternion_priv.h"
#include "pdm_rotation.h"


/*----------------------------------------------------------------------------*/

#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
#ifdef __cplusplus
extern "C"
{
#endif
  void dgemm_(char *transA, char *transB, int *m, int *n, int *k,
                     double *alpha, double *A, int *lda,
                     double *B, int *ldb, double *beta,
                     double *C, int *ldc);
#ifdef __cplusplus
}
#endif
#endif

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

static const double ROTATION_EPS  = __DBL_EPSILON__;

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
      homogeneous_matrix[4*i+j] = (i==j) ? 1. : 0.;
    }
  }
}

static
void
transpose_homogeneous_matrix
(
  double* homogeneous_matrix
)
{
  // transposing the matrix
  //  0  1  2  3
  //  4  5  6  7
  //  8  9 10 11
  // 12 13 14 15
  //
  double tmp;
  int swapped_ind_0[6] = {1,2,6, 3, 7,11};
  int swapped_ind_1[6] = {4,8,9,12,13,14};
  for (int ind=0;ind<6;ind++){
    // swapping values
    tmp = homogeneous_matrix[swapped_ind_0[ind]];
    homogeneous_matrix[swapped_ind_0[ind]] = homogeneous_matrix[swapped_ind_1[ind]];
    homogeneous_matrix[swapped_ind_1[ind]] = tmp;
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
  const PDM_bool_t reverse,
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

// Axis angle to other formats ---

void
PDM_rotation_axis_angle_to_euler_angles
(
  const double axis[3],
  const double angle,
  const PDM_bool_t reverse,
  const int order[3],
  const PDM_bool_t intrinsic,
        double* ang_x,
        double* ang_y,
        double* ang_z
)
{
  _pdm_quaternion_t qt;
  double langle = (reverse ? -angle : angle);
  PDM_quaternion_from_axis_angle(axis,langle,(PDM_quaternion*)&qt);
  PDM_quaternion_to_euler_angles((PDM_quaternion*)&qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void
PDM_rotation_axis_angle_to_rotation_matrix
(
  const double axis[3],
  const double angle,
  const PDM_bool_t reverse,
        double *rotation_matrix
)
{
  _pdm_quaternion_t qt;
  double langle = (reverse ? -angle : angle);
  PDM_quaternion_from_axis_angle(axis,langle,(PDM_quaternion*)&qt);
  PDM_quaternion_to_rotation_matrix((PDM_quaternion*)&qt,rotation_matrix);
}

void
PDM_rotation_axis_angle_to_homogeneous_matrix
(
  const double axis[3],
  const double angle,
  const PDM_bool_t reverse,
        double *homogeneous_matrix
)
{
  _pdm_quaternion_t qt;
  double langle = (reverse ? -angle : angle);
  PDM_quaternion_from_axis_angle(axis,langle,(PDM_quaternion*)&qt);
  PDM_quaternion_to_homogeneous_matrix((PDM_quaternion*)&qt,homogeneous_matrix);
}


void
PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix
(
  const double axis[3],
  const double angle,
  const double rotation_center[3],
  const PDM_bool_t reverse,
        double *homogeneous_matrix
){
  // helper matrices (matrix multiplication (dgemm) is not inplace)
  double trans_mat    [16];
  double rot_mat      [16];
  double rot_trans_mat[16];
  // Trans+.Rot.Trans-
  // translation of -rotation_center
  set_translation_to_homogeneous_matrix(rotation_center,
                                        PDM_TRUE,
                                        trans_mat);
  // rotation
  double langle = (reverse ? -angle : angle);
  PDM_rotation_axis_angle_to_homogeneous_matrix(axis,
                                                langle,
                                                PDM_FALSE,
                                                rot_mat);
  PDM_rotation_multiply_n_by_n_matrices(rot_mat,
                                        trans_mat,
                                        4,
                                        rot_trans_mat);
  // translation of rotation_center
  // re-using trans_mat
  set_translation_to_homogeneous_matrix(rotation_center,
                                        PDM_FALSE,
                                        trans_mat);
  PDM_rotation_multiply_n_by_n_matrices(trans_mat,
                                        rot_trans_mat,
                                        4,
                                        homogeneous_matrix);
}

// Euler angles to other formats ---

void
PDM_rotation_euler_angles_to_axis_angle
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  const PDM_bool_t intrinsic,
  const PDM_bool_t reverse,
        double axis[3],
        double* angle
)
{
  _pdm_quaternion_t qt;
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_axis_angle((PDM_quaternion*)&qt,axis,angle);
}

void
PDM_rotation_euler_angles_to_euler_angles
(
  const double input_ang_x,
  const double input_ang_y,
  const double input_ang_z,
  const int input_order[3],
  const PDM_bool_t input_intrinsic,
  const PDM_bool_t reverse,
  const int output_order[3],
  const PDM_bool_t output_intrinsic,
        double* output_ang_x,
        double* output_ang_y,
        double* output_ang_z
)
{
  _pdm_quaternion_t qt;
  PDM_quaternion_from_euler_angles(input_ang_x,input_ang_y,input_ang_z,input_order,input_intrinsic,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_euler_angles((PDM_quaternion*)&qt,output_order,output_intrinsic,output_ang_x,output_ang_y,output_ang_z);
}

void
PDM_rotation_euler_angles_to_rotation_matrix
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  const PDM_bool_t intrinsic,
  const PDM_bool_t reverse,
        double* rotation_matrix
)
{
  _pdm_quaternion_t qt;
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_rotation_matrix((PDM_quaternion*)&qt,rotation_matrix);
}

void
PDM_rotation_euler_angles_to_homogeneous_matrix
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  const PDM_bool_t intrinsic,
  const PDM_bool_t reverse,
        double* homogeneous_matrix
)
{
  _pdm_quaternion_t qt;
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_homogeneous_matrix((PDM_quaternion*)&qt,homogeneous_matrix);
}

void
PDM_rotation_euler_angles_and_rotation_center_to_homogeneous_matrix
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  const PDM_bool_t intrinsic,
  const double rotation_center[3],
  const PDM_bool_t reverse,
        double* homogeneous_matrix
){
  // helper matrices (matrix multiplication (dgemm) is not inplace)
  double trans_mat    [16];
  double rot_mat      [16];
  double rot_trans_mat[16];
  // Trans+.Rot.Trans-
  // translation of -rotation_center
  set_translation_to_homogeneous_matrix(rotation_center,
                                        PDM_TRUE,
                                        trans_mat);
  // rotation
  PDM_rotation_euler_angles_to_homogeneous_matrix(ang_x,
                                                  ang_y,
                                                  ang_z,
                                                  order,
                                                  intrinsic,
                                                  PDM_FALSE,
                                                  rot_mat);
  if (reverse){
    transpose_homogeneous_matrix(rot_mat);
  }
  PDM_rotation_multiply_n_by_n_matrices(rot_mat,
                                        trans_mat,
                                        4,
                                        rot_trans_mat);
  // translation of rotation_center
  // re-using trans_mat
  set_translation_to_homogeneous_matrix(rotation_center,
                                        PDM_FALSE,
                                        trans_mat);
  PDM_rotation_multiply_n_by_n_matrices(trans_mat,
                                        rot_trans_mat,
                                        4,
                                        homogeneous_matrix);
}

void
PDM_rotation_periodic_t_info_to_homogeneous_matrix
(
  const double rotation_center[3],
  const double rotation_angle[3],
  const double translation[3],
  const PDM_bool_t reverse,
        double *homogeneous_matrix
)
{
  // TODO: raise error if both translation and rotation
  PDM_bool_t apply_translation = (PDM_bool_t) ( PDM_ABS(translation   [0])<ROTATION_EPS || (PDM_ABS(translation   [1])<ROTATION_EPS) || (PDM_ABS(translation   [2])<ROTATION_EPS));
  PDM_bool_t apply_rotation    = (PDM_bool_t) ( PDM_ABS(rotation_angle[0])<ROTATION_EPS || (PDM_ABS(rotation_angle[1])<ROTATION_EPS) || (PDM_ABS(rotation_angle[2])<ROTATION_EPS));
  if (apply_translation & apply_rotation){
    printf("Error : A periodic_t node with both a translation and a rotation is not supported yet.\n");
    exit(EX_DATAERR);
  }
  if (apply_translation){
    // pure translation
    set_identity_to_homogeneous_matrix(homogeneous_matrix);
    set_translation_to_homogeneous_matrix(translation,reverse,homogeneous_matrix);
  }
  else if(apply_rotation) {
    // pure rotation
    const int order[3] = {2,1,0};
    PDM_rotation_euler_angles_and_rotation_center_to_homogeneous_matrix(rotation_angle[0],
                                                                        rotation_angle[1],
                                                                        rotation_angle[2],
                                                                        order,
                                                                        PDM_TRUE,
                                                                        rotation_center,
                                                                        reverse,
                                                                        homogeneous_matrix);
  }
  else {
    set_identity_to_homogeneous_matrix(homogeneous_matrix);
  }
}

// Rotation matrix to other formats ---

void
PDM_rotation_rotation_matrix_to_axis_angle
(
  const double* rotation_matrix,
  const PDM_bool_t reverse,
        double axis[3],
        double* angle
)
{
  _pdm_quaternion_t qt;
  PDM_quaternion_from_rotation_matrix(rotation_matrix,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_axis_angle((PDM_quaternion*)&qt,axis,angle);
}

void
PDM_rotation_rotation_matrix_to_euler_angles
(
  const double* rotation_matrix,
  const PDM_bool_t reverse,
  const int order[3],
  const PDM_bool_t intrinsic,
        double* ang_x,
        double* ang_y,
        double* ang_z
)
{
  _pdm_quaternion_t qt;
  PDM_quaternion_from_rotation_matrix(rotation_matrix,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_euler_angles((PDM_quaternion*)&qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void
PDM_rotation_rotation_matrix_to_homogeneous_matrix
(
  const double* rotation_matrix,
  const PDM_bool_t reverse,
        double* homogeneous_matrix
)
{
  if (reverse){
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        homogeneous_matrix[4*i+j] = rotation_matrix[3*j+i];
      }
    }
  }else{
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        homogeneous_matrix[4*i+j] = rotation_matrix[3*i+j];
      }
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
PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix
(
  const double* rotation_matrix,
  const double  rotation_center[3],
  const PDM_bool_t reverse,
        double* homogeneous_matrix
){
  // helper matrices (matrix multiplication (dgemm) is not inplace)
  double trans_mat    [16];
  double rot_mat      [16];
  double rot_trans_mat[16];
  // Trans+.Rot.Trans-
  // translation of -rotation_center
  set_translation_to_homogeneous_matrix(rotation_center,
                                        PDM_TRUE,
                                        trans_mat);
  // rotation
  PDM_rotation_rotation_matrix_to_homogeneous_matrix(rotation_matrix,
                                                     PDM_FALSE,
                                                     rot_mat);
  if (reverse){
    transpose_homogeneous_matrix(rot_mat);
  }
  PDM_rotation_multiply_n_by_n_matrices(rot_mat,
                                        trans_mat,
                                        4,
                                        rot_trans_mat);
  // translation of rotation_center
  // re-using trans_mat
  set_translation_to_homogeneous_matrix(rotation_center,
                                        PDM_FALSE,
                                        trans_mat);
  PDM_rotation_multiply_n_by_n_matrices(trans_mat,
                                        rot_trans_mat,
                                        4,
                                        homogeneous_matrix);
}

void
PDM_rotation_homogeneous_matrix_to_axis_angle
(
  const double* homogeneous_matrix,
  const PDM_bool_t reverse,
        double axis[3],
        double* angle
)
{
  _pdm_quaternion_t qt;
  PDM_quaternion_from_homogeneous_matrix(homogeneous_matrix,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_axis_angle((PDM_quaternion*)&qt,axis,angle);
}

void
PDM_rotation_homogeneous_matrix_to_euler_angles
(
  const double* homogeneous_matrix,
  const PDM_bool_t reverse,
  const int order[3],
  const PDM_bool_t intrinsic,
        double* ang_x,
        double* ang_y,
        double* ang_z
)
{
  _pdm_quaternion_t qt;
  PDM_quaternion_from_homogeneous_matrix(homogeneous_matrix,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_euler_angles((PDM_quaternion*)&qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void
PDM_rotation_homogeneous_matrix_to_rotation_matrix
(
  const double* homogeneous_matrix,
  const PDM_bool_t reverse,
        double* rotation_matrix
)
{
  if (reverse){
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        rotation_matrix[3*j+i] = homogeneous_matrix[4*i+j];
      }
    }
  }else{
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        rotation_matrix[3*i+j] = homogeneous_matrix[4*i+j];
      }
    }
  }
}

void
PDM_rotation_two_vectors_to_axis_angle
(
  const double vector_1[3],
  const double vector_2[3],
  const PDM_bool_t reverse,
        double axis[3],
        double* angle
)
{
  _pdm_quaternion_t qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_axis_angle((PDM_quaternion*)&qt,axis,angle);
}

void
PDM_rotation_two_vectors_to_euler_angles
(
  const double vector_1[3],
  const double vector_2[3],
  const PDM_bool_t reverse,
  const int order[3],
  const PDM_bool_t intrinsic,
        double* ang_x,
        double* ang_y,
        double* ang_z
)
{
  _pdm_quaternion_t qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_euler_angles((PDM_quaternion*)&qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void
PDM_rotation_two_vectors_to_rotation_matrix
(
  const double vector_1[3],
  const double vector_2[3],
  const PDM_bool_t reverse,
        double *rotation_matrix
)
{
  _pdm_quaternion_t qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_rotation_matrix((PDM_quaternion*)&qt,rotation_matrix);
}

void
PDM_rotation_two_vectors_to_homogeneous_matrix
(
  const double vector_1[3],
  const double vector_2[3],
  const PDM_bool_t reverse,
        double *homogeneous_matrix
)
{
  // Improving precision with quaternion squared
  _pdm_quaternion_t qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,(PDM_quaternion*)&qt);
  if (reverse){
    PDM_quaternion_conjugate((PDM_quaternion*)&qt,(PDM_quaternion*)&qt);
  }
  PDM_quaternion_to_homogeneous_matrix((PDM_quaternion*)&qt,homogeneous_matrix);
  printf("homogeneous_matrix = %12.5e/%12.5e/%12.5e/%12.5e \n", homogeneous_matrix[4*0+0], homogeneous_matrix[4*0+1], homogeneous_matrix[4*0+2], homogeneous_matrix[4*0+3]);
  printf("homogeneous_matrix = %12.5e/%12.5e/%12.5e/%12.5e \n", homogeneous_matrix[4*1+0], homogeneous_matrix[4*1+1], homogeneous_matrix[4*1+2], homogeneous_matrix[4*1+3]);
  printf("homogeneous_matrix = %12.5e/%12.5e/%12.5e/%12.5e \n", homogeneous_matrix[4*2+0], homogeneous_matrix[4*2+1], homogeneous_matrix[4*2+2], homogeneous_matrix[4*2+3]);
  printf("homogeneous_matrix = %12.5e/%12.5e/%12.5e/%12.5e \n", homogeneous_matrix[4*3+0], homogeneous_matrix[4*3+1], homogeneous_matrix[4*3+2], homogeneous_matrix[4*3+3]);
}

void
PDM_rotation_two_vectors_and_rotation_center_to_homogeneous_matrix
(
  const double vector_1[3],
  const double vector_2[3],
  const double rotation_center[3],
  const PDM_bool_t reverse,
        double *homogeneous_matrix
)
{
  // helper matrices (matrix multiplication (dgemm) is not inplace)
  double trans_mat    [16];
  double rot_mat      [16];
  double rot_trans_mat[16];
  // Trans+.Rot.Trans-
  // translation of -rotation_center
  set_translation_to_homogeneous_matrix(rotation_center,
                                        PDM_TRUE,
                                        trans_mat);
  // rotation
  PDM_rotation_two_vectors_to_homogeneous_matrix(vector_1,
                                                 vector_2,
                                                 PDM_FALSE,
                                                 rot_mat);
  if (reverse){
    transpose_homogeneous_matrix(rot_mat);
  }
  PDM_rotation_multiply_n_by_n_matrices(rot_mat,
                                        trans_mat,
                                        4,
                                        rot_trans_mat);
  // translation of rotation_center
  // re-using trans_mat
  set_translation_to_homogeneous_matrix(rotation_center,
                                        PDM_FALSE,
                                        trans_mat);
  PDM_rotation_multiply_n_by_n_matrices(trans_mat,
                                        rot_trans_mat,
                                        4,
                                        homogeneous_matrix);
}

void
PDM_rotation_axes_and_origin_to_homogeneous_matrix
(
  const double axis_1[3],
  const double axis_2[3],
  const double axis_3[3],
  const double origin[3],
  const PDM_bool_t reverse,
        double *homogeneous_matrix
)
{
  double e1_vec_e2[3];
  PDM_bool_t e_is_orthogonal = (PDM_bool_t) ((PDM_ABS(PDM_DOT_PRODUCT(axis_1,axis_2))<ROTATION_EPS) &&\
                                             (PDM_ABS(PDM_DOT_PRODUCT(axis_1,axis_3))<ROTATION_EPS) &&\
                                             (PDM_ABS(PDM_DOT_PRODUCT(axis_2,axis_3))<ROTATION_EPS));
  PDM_CROSS_PRODUCT(e1_vec_e2,axis_1,axis_2);
  PDM_bool_t e_is_direct     = (PDM_bool_t) (PDM_DOT_PRODUCT(e1_vec_e2,axis_3)>0.);
  // normalizing each basis vector
  double inv_norm[3];
  for (int i=0;i<3;i++){
    inv_norm[i] = PDM_DOT_PRODUCT(axis_1,axis_1);
  }
  PDM_bool_t e1_is_null = (PDM_bool_t) (PDM_ABS(inv_norm[0])<ROTATION_EPS);
  PDM_bool_t e2_is_null = (PDM_bool_t) (PDM_ABS(inv_norm[0])<ROTATION_EPS);
  PDM_bool_t e3_is_null = (PDM_bool_t) (PDM_ABS(inv_norm[0])<ROTATION_EPS);
  PDM_bool_t checks_failed = PDM_FALSE;
  if (!e_is_orthogonal){
    printf("Error : Provided axes are not orthogonal: [%10.6e %10.6e %10.6e] // [%10.6e %10.6e %10.6e] // [%10.6e %10.6e %10.6e].\n",axis_1[0],axis_1[1],axis_1[2],axis_2[0],axis_2[1],axis_2[2],axis_3[0],axis_3[1],axis_3[2]);
    checks_failed = PDM_TRUE;
  }
  if (!e_is_direct){
    printf("Error : Provided axes are not direct.\n");
    checks_failed = PDM_TRUE;
  }
  if (e1_is_null){
    printf("Error : First axis vector [%10.6e %10.6e %10.6e] is null (norm squared = %.16e).\n",axis_1[0],axis_1[1],axis_1[2],inv_norm[0]);
    checks_failed = PDM_TRUE;
  }
  if (e2_is_null){
    printf("Error : Second axis vector [%10.6e %10.6e %10.6e] is null (norm squared = %.16e).\n",axis_2[0],axis_2[1],axis_2[2],inv_norm[1]);
    checks_failed = PDM_TRUE;
  }
  if (e3_is_null){
    printf("Error : Third axis vector [%10.6e %10.6e %10.6e] is null (norm squared = %.16e).\n",axis_3[0],axis_3[1],axis_3[2],inv_norm[2]);
    checks_failed = PDM_TRUE;
  }
  if (checks_failed){
    exit(EX_USAGE);
  }
  for (int i=0;i<3;i++){
    inv_norm[i] = 1./sqrt(inv_norm[i]);
  }
  double translation_vec[3];
  set_identity_to_homogeneous_matrix(homogeneous_matrix);
  if (reverse){
    // with homogeneous matrices the translation occurs in the second reference frame:
    // y = H.x = R.x + V
    // x = R^T.y + V_rev
    // y = R.(R^T.y + V_rev) + V = y + R.V_rev + V
    // thus R.V_rev = - V => V_rev = -R^T.V
    // building R^T
    double rot_mat_T[9];
    rot_mat_T[3*0+0] = inv_norm[0]*axis_1[0];
    rot_mat_T[3*0+1] = inv_norm[0]*axis_1[1];
    rot_mat_T[3*0+2] = inv_norm[0]*axis_1[2];
    rot_mat_T[3*1+0] = inv_norm[1]*axis_2[0];
    rot_mat_T[3*1+1] = inv_norm[1]*axis_2[1];
    rot_mat_T[3*1+2] = inv_norm[1]*axis_2[2];
    rot_mat_T[3*2+0] = inv_norm[2]*axis_3[0];
    rot_mat_T[3*2+1] = inv_norm[2]*axis_3[1];
    rot_mat_T[3*2+2] = inv_norm[2]*axis_3[2];
    PDM_rotation_apply_n_by_n_matrix(rot_mat_T,origin,3,1,translation_vec);
    homogeneous_matrix[4*0+3] = -translation_vec[0];
    homogeneous_matrix[4*1+3] = -translation_vec[1];
    homogeneous_matrix[4*2+3] = -translation_vec[2];
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        homogeneous_matrix[4*i+j] = rot_mat_T[3*i+j];
      }
    }
  }else{
    homogeneous_matrix[4*0+0] = inv_norm[0]*axis_1[0];
    homogeneous_matrix[4*0+1] = inv_norm[1]*axis_2[0];
    homogeneous_matrix[4*0+2] = inv_norm[2]*axis_3[0];
    homogeneous_matrix[4*1+0] = inv_norm[0]*axis_1[1];
    homogeneous_matrix[4*1+1] = inv_norm[1]*axis_2[1];
    homogeneous_matrix[4*1+2] = inv_norm[2]*axis_3[1];
    homogeneous_matrix[4*2+0] = inv_norm[0]*axis_1[2];
    homogeneous_matrix[4*2+1] = inv_norm[1]*axis_2[2];
    homogeneous_matrix[4*2+2] = inv_norm[2]*axis_3[2];
    homogeneous_matrix[4*0+3] = origin[0];
    homogeneous_matrix[4*1+3] = origin[1];
    homogeneous_matrix[4*2+3] = origin[2];
  }
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
  char no_trans = 'N';
  double alpha = 1.;
  double beta = 0.;
  // avoiding const warnings
  int n_ = n;
  double* A_ = (double*) A;
  double* B_ = (double*) B;
  dgemm_(&no_trans,
         &no_trans,
         &n_,
         &n_,
         &n_,
         &alpha,
         B_,
         &n_,
         A_,
         &n_,
         &beta,
         C,
         &n_);
#else
  PDM_UNUSED(A);
  PDM_UNUSED(B);
  PDM_UNUSED(n);
  PDM_UNUSED(C);
  printf("Error : BLAS is mandatory (shipped with LAPACK or MKL), recompile with it.\n");
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
  // avoiding const warnings
  char no_trans = 'N';
  char    trans = 'T';
  double alpha = 1.;
  double beta  = 0.;
  double* A_ = (double*) A; // (n-by-n) -> k=n / n=n
  double* x_ = (double*) x; // (n_samp-by-n) -> m=n_samp / k=n
  int m_ = n_samp;
  int n_ = n;
  int k_ = n; // sure
  dgemm_(&trans,
         &no_trans,
         &n_,       // n (flipping m & n > fortran/C)
         &m_,       // m
         &k_,       // k
         &alpha,
         A_,
         &n_,       // lda
         x_,
         &n_,       // ldb
         &beta,
         y,
         &n_);      // ldc
#else
  PDM_UNUSED(A);
  PDM_UNUSED(x);
  PDM_UNUSED(n);
  PDM_UNUSED(n_samp);
  PDM_UNUSED(y);
  printf("Error : BLAS is mandatory (shipped with LAPACK or MKL), recompile with it.\n");
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
  double homogeneous_matrix    [16];
  PDM_rotation_euler_angles_and_rotation_center_to_homogeneous_matrix(ang_x,
                                                                      ang_y,
                                                                      ang_z,
                                                                      order,
                                                                      intrinsic,
                                                                      rotation_center,
                                                                      reverse,
                                                                      homogeneous_matrix);
  // // building the homogeneous matrix corresponding to whole transformation:
  // // Trans+.Rot.Trans-
  // double homogeneous_matrix_tmp[16];
  // double tmp_matrix[16];
  // set_identity_to_homogeneous_matrix(homogeneous_matrix);
  // // translation of -rotation_center
  // set_translation_to_homogeneous_matrix(rotation_center,PDM_TRUE,
  //   homogeneous_matrix);

  // // rotation
  // PDM_rotation_euler_angles_to_homogeneous_matrix(ang_x,ang_y,ang_z,order,
  //   intrinsic,tmp_matrix);
  // if (reverse){
  //   transpose_homogeneous_matrix(tmp_matrix);
  // }
  // PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
  //                                       homogeneous_matrix,
  //                                       4,
  //                                       homogeneous_matrix_tmp);
  // for(int i = 0; i < 16; ++i) {
  //   homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  // }

  // // translation of rotation_center
  // set_translation_to_homogeneous_matrix(rotation_center,
  //                                       PDM_FALSE,
  //                                       tmp_matrix);
  // PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
  //                                       homogeneous_matrix,
  //                                       4,
  //                                       homogeneous_matrix_tmp);
  // for(int i = 0; i < 16; ++i) {
  //   homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  // }

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
  double homogeneous_matrix    [16];
  PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix(axis,
                                                                    angle,
                                                                    rotation_center,
                                                                    reverse,
                                                                    homogeneous_matrix);
  // // building the homogeneous matrix corresponding to whole transformation:
  // // Trans+.Rot.Trans-
  // double homogeneous_matrix[16];
  // double homogeneous_matrix_tmp[16];
  // double tmp_matrix[16];
  // set_identity_to_homogeneous_matrix(homogeneous_matrix);
  // // translation of -rotation_center
  // set_translation_to_homogeneous_matrix(rotation_center,PDM_TRUE,
  //   homogeneous_matrix);

  // // rotation
  // double langle = (reverse ? -angle : angle);
  // PDM_rotation_axis_angle_to_homogeneous_matrix(axis,langle,tmp_matrix);

  // // BLAS_DGEMM does not support inplace
  // PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
  //                                       homogeneous_matrix,
  //                                       4,
  //                                       homogeneous_matrix_tmp);

  // for(int i = 0; i < 16; ++i) {
  //   homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  // }

  // // translation of rotation_center
  // set_translation_to_homogeneous_matrix(rotation_center,
  //                                       PDM_FALSE,
  //                                       tmp_matrix);

  // PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
  //                                       homogeneous_matrix,
  //                                       4,
  //                                       homogeneous_matrix_tmp);

  // for(int i = 0; i < 16; ++i) {
  //   homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  // }

  // applying the homogeneous matrix to the coordinate vector
  PDM_rotation_apply_homogeneous_matrix(homogeneous_matrix,
                                        vector,
                                        n_samp,
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
  double homogeneous_matrix[16];
  PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix(rotation_matrix,
                                                                         rotation_center,
                                                                         reverse,
                                                                         homogeneous_matrix);
  // // building the homogeneous matrix corresponding to whole transformation:
  // // Trans+.Rot.Trans-
  // double homogeneous_matrix_tmp[16];
  // double tmp_matrix[16];
  // // translation of -rotation_center
  // set_translation_to_homogeneous_matrix(rotation_center,PDM_TRUE,
  //   homogeneous_matrix);

  // // rotation
  // PDM_rotation_rotation_matrix_to_homogeneous_matrix(rotation_matrix,
  //   tmp_matrix);
  // if (reverse){
  //   transpose_homogeneous_matrix(tmp_matrix);
  // }
  // PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
  //                                       homogeneous_matrix,
  //                                       4,
  //                                       homogeneous_matrix_tmp);

  // for(int i = 0; i < 16; ++i) {
  //   homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  // }

  // // translation of rotation_center
  // set_translation_to_homogeneous_matrix(rotation_center,
  //                                       PDM_FALSE,
  //                                       tmp_matrix);

  // PDM_rotation_multiply_n_by_n_matrices(tmp_matrix,
  //                                       homogeneous_matrix,
  //                                       4,
  //                                       homogeneous_matrix_tmp);

  // for(int i = 0; i < 16; ++i) {
  //   homogeneous_matrix[i] = homogeneous_matrix_tmp[i];
  // }

  // applying the homogeneous matrix to the coordinate vector
  PDM_rotation_apply_homogeneous_matrix(homogeneous_matrix,vector,n_samp,
    vector_out);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
