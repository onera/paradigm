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
void
dgemm_
(
  char   *transA,
  char   *transB,
  int    *m,
  int    *n,
  int    *k,
  double *alpha,
  double *A,
  int    *lda,
  double *B,
  int    *ldb,
  double *beta,
  double *C,
  int    *ldc
);
#endif


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
  PDM_quaternion qt;
  double langle = (reverse ? -angle : angle);
  PDM_quaternion_from_axis_angle(axis,langle,&qt);
  PDM_quaternion_to_euler_angles(&qt,order,intrinsic,ang_x,ang_y,ang_z);
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
  PDM_quaternion qt;
  double langle = (reverse ? -angle : angle);
  PDM_quaternion_from_axis_angle(axis,langle,&qt);
  PDM_quaternion_to_rotation_matrix(&qt,rotation_matrix);
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
  PDM_quaternion qt;
  double langle = (reverse ? -angle : angle);
  PDM_quaternion_from_axis_angle(axis,langle,&qt);
  PDM_quaternion_to_homogeneous_matrix(&qt,homogeneous_matrix);
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
  double rot_mat[9];
  // rotation
  double langle = (reverse ? -angle : angle);
  PDM_rotation_axis_angle_to_rotation_matrix(axis,
                                             langle,
                                             PDM_FALSE,
                                             rot_mat);
  PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix(rot_mat,
                                                                         rotation_center,
                                                                         PDM_FALSE,
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
  PDM_quaternion qt;
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_axis_angle(&qt,axis,angle);
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
  PDM_quaternion qt;
  PDM_quaternion_from_euler_angles(input_ang_x,input_ang_y,input_ang_z,input_order,input_intrinsic,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_euler_angles(&qt,output_order,output_intrinsic,output_ang_x,output_ang_y,output_ang_z);
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
  PDM_quaternion qt;
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_rotation_matrix(&qt,rotation_matrix);
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
  PDM_quaternion qt;
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_homogeneous_matrix(&qt,homogeneous_matrix);
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
  PDM_bool_t apply_translation = (PDM_bool_t) ( PDM_ABS(translation   [0])>ROTATION_EPS || (PDM_ABS(translation   [1])>ROTATION_EPS) || (PDM_ABS(translation   [2])>ROTATION_EPS));
  PDM_bool_t apply_rotation    = (PDM_bool_t) ( PDM_ABS(rotation_angle[0])>ROTATION_EPS || (PDM_ABS(rotation_angle[1])>ROTATION_EPS) || (PDM_ABS(rotation_angle[2])>ROTATION_EPS));
  if (apply_translation & apply_rotation){
    const int order[3] = {2,1,0};
    // forcing rotation center to 0, rotation center is considered below with translation
    PDM_rotation_euler_angles_to_homogeneous_matrix(rotation_angle[0],
                                                    rotation_angle[1],
                                                    rotation_angle[2],
                                                    order,
                                                    PDM_TRUE,
                                                    reverse,
                                                    homogeneous_matrix);
    // if C is the rotation center, T the provided translation and R the rotation matrix
    // homogenous translation is tau:
    double RC[3]; // R.C
    PDM_rotation_apply_homogeneous_matrix(homogeneous_matrix,rotation_center,1,RC);
    if (reverse){
      // tau = C - R.C - R.T (reverse) (note that R is computed as R^-1)
      double RT[3]; // R.T
      PDM_rotation_apply_homogeneous_matrix(homogeneous_matrix,translation,1,RT);
      homogeneous_matrix[3]  = rotation_center[0] - RC[0] - RT[0];
      homogeneous_matrix[7]  = rotation_center[1] - RC[1] - RT[1];
      homogeneous_matrix[11] = rotation_center[2] - RC[2] - RT[2];
    }
    else{
      // tau = C - R.C + T (direct)
      homogeneous_matrix[3]  = rotation_center[0] - RC[0] + translation[0];
      homogeneous_matrix[7]  = rotation_center[1] - RC[1] + translation[1];
      homogeneous_matrix[11] = rotation_center[2] - RC[2] + translation[2];
    }
  }
  else if (apply_translation){
    // pure translation
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
  PDM_quaternion qt;
  PDM_quaternion_from_rotation_matrix(rotation_matrix,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_axis_angle(&qt,axis,angle);
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
  PDM_quaternion qt;
  PDM_quaternion_from_rotation_matrix(rotation_matrix,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_euler_angles(&qt,order,intrinsic,ang_x,ang_y,ang_z);
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
                                                     reverse,
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

void
PDM_rotation_homogeneous_matrix_to_axis_angle
(
  const double* homogeneous_matrix,
  const PDM_bool_t reverse,
        double axis[3],
        double* angle
)
{
  PDM_quaternion qt;
  PDM_quaternion_from_homogeneous_matrix(homogeneous_matrix,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_axis_angle(&qt,axis,angle);
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
  PDM_quaternion qt;
  PDM_quaternion_from_homogeneous_matrix(homogeneous_matrix,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_euler_angles(&qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void
PDM_rotation_homogeneous_matrix_to_euler_angles_and_translation
(
  const double* homogeneous_matrix,
  const PDM_bool_t reverse,
  const int order[3],
  const PDM_bool_t intrinsic,
  double* ang_x,
  double* ang_y,
  double* ang_z,
  double  translation[3]
)
{
  PDM_quaternion qt;
  PDM_quaternion_from_homogeneous_matrix(homogeneous_matrix,&qt);
  double axis[3];
  double angle;
  PDM_quaternion_to_axis_angle(&qt,axis,&angle);
  if (PDM_ABS(angle)<ROTATION_EPS){
    // no rotation
    (*ang_x) = 0.;
    (*ang_y) = 0.;
    (*ang_z) = 0.;
    if (reverse){
      translation[0] = -homogeneous_matrix[3];
      translation[1] = -homogeneous_matrix[7];
      translation[2] = -homogeneous_matrix[11];
    }else{
      translation[0] = homogeneous_matrix[3];
      translation[1] = homogeneous_matrix[7];
      translation[2] = homogeneous_matrix[11];
    }
  } else {
    // Y = RX + T (direct)
    // X = RY + T (inverse) -> Y = R^T.X - R^T.T
    PDM_rotation_homogeneous_matrix_to_euler_angles(homogeneous_matrix,reverse,order,intrinsic,ang_x,ang_y,ang_z);
    // // Translation part of homogeneous matrix T is:
    // // T = C - R.C + T' (C rotation center, R rotation matrix, T' remaining translation)
    // // C - R.C = T - T' = T - (T.axis)axis
    if (reverse){
      double rot_mat_inv[9];
      PDM_rotation_homogeneous_matrix_to_rotation_matrix(homogeneous_matrix,PDM_TRUE,rot_mat_inv);
      double T[3] = {
        homogeneous_matrix[3],
        homogeneous_matrix[7],
        homogeneous_matrix[11]
      };
      PDM_rotation_apply_n_by_n_matrix(rot_mat_inv,T,3,1,translation);
      translation[0] *= -1.;
      translation[1] *= -1.;
      translation[2] *= -1.;
    }else{
      translation[0] = homogeneous_matrix[3];
      translation[1] = homogeneous_matrix[7];
      translation[2] = homogeneous_matrix[11];
    }
  }
}

void
PDM_rotation_homogeneous_matrix_to_periodic_t_info
(
  const double* homogeneous_matrix,
  const PDM_bool_t reverse,
  const PDM_bool_t compute_rotation_center,
        double rotation_center[3],
        double rotation_angle[3],
        double translation[3]
)
{
  const int order[3] = {2,1,0};
  double ang_x,ang_y,ang_z;
  if (compute_rotation_center){
    PDM_quaternion qt;
    PDM_quaternion_from_homogeneous_matrix(homogeneous_matrix,&qt);
    double axis[3];
    double angle;
    double sign = 1.;
    if (reverse)
      sign = -1.;
    PDM_quaternion_to_axis_angle(&qt,axis,&angle);
    if (PDM_ABS(angle)<ROTATION_EPS){
      // no rotation
      ang_x = 0.;
      ang_y = 0.;
      ang_z = 0.;
      translation[0]     = homogeneous_matrix[3]*sign;
      translation[1]     = homogeneous_matrix[7]*sign;
      translation[2]     = homogeneous_matrix[11]*sign;
      rotation_center[0] = 0.;
      rotation_center[1] = 0.;
      rotation_center[2] = 0.;

    } else {
      PDM_rotation_homogeneous_matrix_to_euler_angles(homogeneous_matrix,
                                                      reverse,
                                                      order,
                                                      PDM_TRUE,
                                                      &ang_x,
                                                      &ang_y,
                                                      &ang_z);

      // Compute change in reference frame (x',y',z') to have z' aligned with the rotation axis
      const double new_z[3] = {0.,0.,1.};
      PDM_quaternion qt_axis;
      PDM_quaternion_from_two_vectors(axis,new_z,&qt_axis);
      // now working in (x',y') of the new reference frame
      double mat_direct[16];
      double mat_inverse[16];
      PDM_quaternion_to_homogeneous_matrix(&qt_axis,mat_direct);
      PDM_quaternion_to_homogeneous_matrix(&qt_axis,mat_inverse);
      transpose_homogeneous_matrix(mat_inverse);
      const double* matrix_list[3] = {mat_direct,(double*)homogeneous_matrix,mat_inverse};
      double homo_mat_in_new_frame[16];
      PDM_rotation_compose_homogeneous_matrices(matrix_list,3,homo_mat_in_new_frame);
      // homo_mat_in_new_frame has now non-identity components on (x',y') only
      // in the reference frame, affine components on x' & y' are:
      // C-R.C = (Id-R).C -> solving for C
      double C_minus_RC[2] = {homo_mat_in_new_frame[3],homo_mat_in_new_frame[7]};
      double id_minus_R_adj[4] = { // adjoint of Id-R
        1.-homo_mat_in_new_frame[0],    homo_mat_in_new_frame[1],
           homo_mat_in_new_frame[4], 1.-homo_mat_in_new_frame[5],
      };
      double inv_det_id_minus_R_adj = 1./(id_minus_R_adj[0]*id_minus_R_adj[3] - id_minus_R_adj[1]*id_minus_R_adj[2]);
      double rotation_center_in_new_frame[3] = {
        inv_det_id_minus_R_adj*(id_minus_R_adj[0]*C_minus_RC[0]+id_minus_R_adj[1]*C_minus_RC[1]),
        inv_det_id_minus_R_adj*(id_minus_R_adj[2]*C_minus_RC[0]+id_minus_R_adj[3]*C_minus_RC[1]),
        0.
      };
      // reverting to old reference frame (x,y,z)
      PDM_rotation_apply_homogeneous_matrix(mat_inverse,rotation_center_in_new_frame,1,rotation_center);
      // Computing translation T' as:
      double T[3] = {
        homogeneous_matrix[3],
        homogeneous_matrix[7],
        homogeneous_matrix[11]
      };
      if (reverse){
        // T' = C - R.C - T
        double R[9];
        double RC[3];
        PDM_rotation_homogeneous_matrix_to_rotation_matrix(homogeneous_matrix,PDM_FALSE,R);
        PDM_rotation_apply_n_by_n_matrix(R,rotation_center,3,1,RC);
        translation[0] = rotation_center[0] - RC[0] - T[0];
        translation[1] = rotation_center[1] - RC[1] - T[1];
        translation[2] = rotation_center[2] - RC[2] - T[2];
      }else{
        // T' = T -C + R.C = T + (R-Id).C (direct)
        double R_minus_id[9];
        PDM_rotation_homogeneous_matrix_to_rotation_matrix(homogeneous_matrix,PDM_FALSE,R_minus_id);
        R_minus_id[0] -= 1.;
        R_minus_id[4] -= 1.;
        R_minus_id[8] -= 1.;
        PDM_rotation_apply_n_by_n_matrix(R_minus_id,rotation_center,3,1,translation); // (R-Id).C
        translation[0] += T[0];
        translation[1] += T[1];
        translation[2] += T[2];
      }
    }
  } else {
    PDM_rotation_homogeneous_matrix_to_euler_angles_and_translation(homogeneous_matrix,
                                                                    reverse,
                                                                    order,
                                                                    PDM_TRUE,
                                                                    &ang_x,
                                                                    &ang_y,
                                                                    &ang_z,
                                                                    translation);
    rotation_center[0] = 0.;
    rotation_center[1] = 0.;
    rotation_center[2] = 0.;
  }
  rotation_angle[0]  = ang_x;
  rotation_angle[1]  = ang_y;
  rotation_angle[2]  = ang_z;
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
  PDM_quaternion qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_axis_angle(&qt,axis,angle);
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
  PDM_quaternion qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_euler_angles(&qt,order,intrinsic,ang_x,ang_y,ang_z);
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
  PDM_quaternion qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_rotation_matrix(&qt,rotation_matrix);
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
  PDM_quaternion qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,&qt);
  if (reverse){
    PDM_quaternion_conjugate(&qt,&qt);
  }
  PDM_quaternion_to_homogeneous_matrix(&qt,homogeneous_matrix);
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
  double rot_mat[9];
  PDM_rotation_two_vectors_to_rotation_matrix(vector_1,
                                              vector_2,
                                              reverse,
                                              rot_mat);
  PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix(rot_mat,
                                                                         rotation_center,
                                                                         PDM_FALSE,
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
  }else{
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
    homogeneous_matrix[4*3+3] = 1.;
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        homogeneous_matrix[4*i+j] = rot_mat_T[3*i+j];
      }
    }
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
  double rotation_matrix[9];
  PDM_rotation_homogeneous_matrix_to_rotation_matrix(homogeneous_matrix,PDM_FALSE,rotation_matrix);
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
  double output_matrix_tmp[16];
  set_identity_to_homogeneous_matrix(output_matrix);
  for (int i = 0; i < n_matrices; i++) {
    PDM_rotation_multiply_n_by_n_matrices(output_matrix,
                                          homogeneous_matrices[i],
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

  // applying the homogeneous matrix to the coordinate vector
  PDM_rotation_apply_homogeneous_matrix(homogeneous_matrix,vector,n_samp,
    vector_out);
}

