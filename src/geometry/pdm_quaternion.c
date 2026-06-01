/*============================================================================
 * Mesure des temps CPU et elapsed
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_quaternion.h"
#include "pdm_quaternion_priv.h"

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

static const double QUATERNION_EPS  = __DBL_EPSILON__;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_PDM_quaternion_compute_squared
(
  PDM_quaternion* qt
)
{
  qt->q_squared[0] = qt->q[0] * qt->q[0];
  qt->q_squared[1] = qt->q[1] * qt->q[1];
  qt->q_squared[2] = qt->q[2] * qt->q[2];
  qt->q_squared[3] = qt->q[3] * qt->q[3];
}

// static
// double
// _PDM_quaternion_3x3_mat_det
// (
//   const double rot_mat[9]
// )
// { return (rot_mat[3*0+0] * (rot_mat[3*1+1] * rot_mat[3*2+2] - rot_mat[3*2+1] * rot_mat[3*1+2])
//          -rot_mat[3*1+0] * (rot_mat[3*0+1] * rot_mat[3*2+2] - rot_mat[3*2+1] * rot_mat[3*0+2])
//          +rot_mat[3*2+0] * (rot_mat[3*0+1] * rot_mat[3*1+2] - rot_mat[3*1+1] * rot_mat[3*0+2]));
// }

/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_quaternion_set
(
  const double w,
  const double v1,
  const double v2,
  const double v3,
        PDM_quaternion* qt
)
{
  qt->q[0] = w;
  qt->q[1] = v1;
  qt->q[2] = v2;
  qt->q[3] = v3;
  _PDM_quaternion_compute_squared(qt);
}


void
PDM_quaternion_set_identity
(
  PDM_quaternion* qt
)
{
  qt->q[0] = 1.;
  qt->q[1] = 0.;
  qt->q[2] = 0.;
  qt->q[3] = 0.;
  qt->q_squared[0] = 1.;
  qt->q_squared[1] = 0.;
  qt->q_squared[2] = 0.;
  qt->q_squared[3] = 0.;
}

void
PDM_quaternion_print
(
  const PDM_quaternion* qt
)
{
  // PDM_printf("qt = [%23.16e,%23.16e,%23.16e,%23.16e]\n",
  PDM_printf("qt = [%12.5e %12.5e %12.5e %12.5e]\n",
    qt->q[0],
    qt->q[1],
    qt->q[2],
    qt->q[3]);
  }

  void
  PDM_quaternion_conjugate
  (
    const PDM_quaternion* qt,
    PDM_quaternion* qt_out
  )
  {
  qt_out->q[0] =  qt->q[0];
  qt_out->q[1] = -qt->q[1];
  qt_out->q[2] = -qt->q[2];
  qt_out->q[3] = -qt->q[3];
  qt_out->q_squared[0] = qt->q_squared[0];
  qt_out->q_squared[1] = qt->q_squared[1];
  qt_out->q_squared[2] = qt->q_squared[2];
  qt_out->q_squared[3] = qt->q_squared[3];
}

double
PDM_quaternion_norm
(
  const PDM_quaternion* qt
)
{
  return sqrt(qt->q_squared[0]+qt->q_squared[1]+qt->q_squared[2]+qt->q_squared[3]);
}

void
PDM_quaternion_normalize
(
  PDM_quaternion* qt
)
{
  double inv_len    = 1. / PDM_quaternion_norm(qt);
  double inv_len_sq = 1. / (qt->q_squared[0]+qt->q_squared[1]+qt->q_squared[2]+qt->q_squared[3]);
  qt->q[0] *= inv_len;
  qt->q[1] *= inv_len;
  qt->q[2] *= inv_len;
  qt->q[3] *= inv_len;
  qt->q_squared[0] *= inv_len_sq;
  qt->q_squared[1] *= inv_len_sq;
  qt->q_squared[2] *= inv_len_sq;
  qt->q_squared[3] *= inv_len_sq;
}

void
PDM_quaternion_compose
(
  const PDM_quaternion* qt_1,
  const PDM_quaternion* qt_2,
  PDM_quaternion* qt_out
)
{
  // buffering data in case where qt_1==qt_out or qt_2==qt_out
  double w,v1,v2,v3;
  w  = qt_1->q[0] * qt_2->q[0] - qt_1->q[1] * qt_2->q[1] - qt_1->q[2] * qt_2->q[2] - qt_1->q[3] * qt_2->q[3];
  v1 = qt_1->q[1] * qt_2->q[0] + qt_1->q[0] * qt_2->q[1] + qt_1->q[2] * qt_2->q[3] - qt_1->q[3] * qt_2->q[2];
  v2 = qt_1->q[0] * qt_2->q[2] - qt_1->q[1] * qt_2->q[3] + qt_1->q[2] * qt_2->q[0] + qt_1->q[3] * qt_2->q[1];
  v3 = qt_1->q[0] * qt_2->q[3] + qt_1->q[1] * qt_2->q[2] - qt_1->q[2] * qt_2->q[1] + qt_1->q[3] * qt_2->q[0] ;
  qt_out->q[0] = w ;
  qt_out->q[1] = v1;
  qt_out->q[2] = v2;
  qt_out->q[3] = v3;
  _PDM_quaternion_compute_squared(qt_out);
  // _qt_out->q_squared[0] = w  * w;
  // _qt_out->q_squared[1] = v1 * v1;
  // _qt_out->q_squared[2] = v2 * v2;
  // _qt_out->q_squared[3] = v3 * v3;
}

PDM_bool_t
PDM_quaternion_equal
(
  const PDM_quaternion* qt,
  const double w,
  const double v0,
  const double v1,
  const double v2,
  const double epsilon
)
{

  return (PDM_bool_t) ((PDM_ABS(qt->q[0] - w ) < epsilon) && \
  (PDM_ABS(qt->q[1] - v0) < epsilon) && \
  (PDM_ABS(qt->q[2] - v1) < epsilon) && \
  (PDM_ABS(qt->q[3] - v2) < epsilon));
}

PDM_bool_t
PDM_quaternion_equal_quaternion
(
  const PDM_quaternion* qt_1,
  const PDM_quaternion* qt_2,
        double          epsilon
)
{
  return PDM_quaternion_equal(qt_1,qt_2->q[0],qt_2->q[1],qt_2->q[2],qt_2->q[3],epsilon);
}


/*----------------------------------------------------------------------------
*  APPLICATION TO VECTORS
*----------------------------------------------------------------------------*/

void
PDM_quaternion_rotate
(
  const PDM_quaternion* qt,
  const double*         vector,
  const int             n_samp,
  double*         vector_out
)
{
  double ww = qt->q_squared[0];
  double xx = qt->q_squared[1];
  double yy = qt->q_squared[2];
  double zz = qt->q_squared[3];
  double wx = qt->q[0] * qt->q[1];
  double wy = qt->q[0] * qt->q[2];
  double wz = qt->q[0] * qt->q[3];
  double xy = qt->q[1] * qt->q[2];
  double xz = qt->q[1] * qt->q[3];
  double yz = qt->q[2] * qt->q[3];

  for (int i = 0; i < n_samp; i++) {
    vector_out[3*i+0] =   ww*vector[3*i+0] + 2*wy*vector[3*i+2] - 2*wz*vector[3*i+1] +   xx*vector[3*i+0] +\
                        2*xy*vector[3*i+1] + 2*xz*vector[3*i+2] -   zz*vector[3*i+0] -   yy*vector[3*i+0];
    vector_out[3*i+1] = 2*xy*vector[3*i+0] +   yy*vector[3*i+1] + 2*yz*vector[3*i+2] + 2*wz*vector[3*i+0] -\
                          zz*vector[3*i+1] +   ww*vector[3*i+1] - 2*wx*vector[3*i+2] -   xx*vector[3*i+1];
    vector_out[3*i+2] = 2*xz*vector[3*i+0] + 2*yz*vector[3*i+1] +   zz*vector[3*i+2] - 2*wy*vector[3*i+0] -\
                          yy*vector[3*i+2] + 2*wx*vector[3*i+1] -   xx*vector[3*i+2] +   ww*vector[3*i+2];
  }
}

void
PDM_quaternion_rotate_derivative
(
  const PDM_quaternion* qt,
  const PDM_quaternion* qt_der,
  const double* vector,
  const int n_samp,
  double* vector_output_der
)
{
  // /!\ not tested in pdm_quaternion.test.cpp yet !
  double ww_der = 2*qt->q[0] * qt_der->q[0];
  double xx_der = 2*qt->q[1] * qt_der->q[1];
  double yy_der = 2*qt->q[2] * qt_der->q[2];
  double zz_der = 2*qt->q[3] * qt_der->q[3];
  double wx_der = qt_der->q[0] * qt->q[1] + qt->q[0] * qt_der->q[1];
  double wy_der = qt_der->q[0] * qt->q[2] + qt->q[0] * qt_der->q[2];
  double wz_der = qt_der->q[0] * qt->q[3] + qt->q[0] * qt_der->q[3];
  double xy_der = qt_der->q[1] * qt->q[2] + qt->q[1] * qt_der->q[2];
  double xz_der = qt_der->q[1] * qt->q[3] + qt->q[1] * qt_der->q[3];
  double yz_der = qt_der->q[2] * qt->q[3] + qt->q[2] * qt_der->q[3];

  for (int i = 0; i < n_samp; i++) {
    vector_output_der[3*i+0] =   ww_der*vector[3*i+0] + 2*wy_der*vector[3*i+2] - 2*wz_der*vector[3*i+1] +   xx_der*vector[3*i+0] +\
                               2*xy_der*vector[3*i+1] + 2*xz_der*vector[3*i+2] -   zz_der*vector[3*i+0] -   yy_der*vector[3*i+0];
    vector_output_der[3*i+1] = 2*xy_der*vector[3*i+0] +   yy_der*vector[3*i+1] + 2*yz_der*vector[3*i+2] + 2*wz_der*vector[3*i+0] -\
                                 zz_der*vector[3*i+1] +   ww_der*vector[3*i+1] - 2*wx_der*vector[3*i+2] -   xx_der*vector[3*i+1];
    vector_output_der[3*i+2] = 2*xz_der*vector[3*i+0] + 2*yz_der*vector[3*i+1] +   zz_der*vector[3*i+2] - 2*wy_der*vector[3*i+0] -\
                                 yy_der*vector[3*i+2] + 2*wx_der*vector[3*i+1] -   xx_der*vector[3*i+2] +   ww_der*vector[3*i+2];
  }
}

/*----------------------------------------------------------------------------
 *  SLERP
 *----------------------------------------------------------------------------*/

void
PDM_quaternion_slerp_from_two_vectors
(
  const double vector_1[3],
  const double vector_2[3],
  const double* t,
  const int n_samp,
        double* out
)
{
  // /!\ not tested in pdm_quaternion.test.cpp yet !
  PDM_quaternion qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,&qt);
  double axis[3];
  double angle;
  PDM_quaternion_to_axis_angle(&qt,axis,&angle);
  for (int i = 0; i < n_samp; i++) {
    PDM_quaternion_from_axis_angle(axis,t[i]*angle,&qt);
    PDM_quaternion_rotate(&qt,vector_1,1,&out[3*i]);
  }
}

void
PDM_quaternion_slerp_from_two_vectors_derivative
(
  const double  vector_1[3],
  const double  vector_2[3],
  const double* t,
  const int     n_samp,
        double* out
)
{
  // /!\ not tested in pdm_quaternion.test.cpp yet !
  PDM_quaternion qt;
  PDM_quaternion qt_der;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,&qt);
  double axis[3];
  double axis_der[3] = {0,0,0};
  double angle;
  PDM_quaternion_to_axis_angle(&qt,axis,&angle);
  for (int i = 0; i < n_samp; i++) {
    PDM_quaternion_from_axis_angle(axis,t[i]*angle,&qt);
    PDM_quaternion_from_axis_angle_derivative(axis,t[i]*angle,axis_der,angle,&qt_der);
    PDM_quaternion_rotate_derivative(&qt,&qt_der,vector_1,1,&out[3*i]);
  }
}


/*----------------------------------------------------------------------------
*  PDM_quaternion_t FROM OTHER FORMATS
*----------------------------------------------------------------------------*/

// https://stackoverflow.com/questions/1171849/finding-quaternion-representing-the-rotation-from-one-vector-to-another
void
PDM_quaternion_from_two_vectors
(
  const double*         vector_1,
  const double*         vector_2,
        PDM_quaternion* qt_out
)
{
  double normsq_1 = vector_1[0]*vector_1[0] + vector_1[1]*vector_1[1] + vector_1[2]*vector_1[2];
  double normsq_2 = vector_2[0]*vector_2[0] + vector_2[1]*vector_2[1] + vector_2[2]*vector_2[2];
  double norm_1 = sqrt(vector_1[0]*vector_1[0] + vector_1[1]*vector_1[1] + vector_1[2]*vector_1[2]);
  double norm_2 = sqrt(vector_2[0]*vector_2[0] + vector_2[1]*vector_2[1] + vector_2[2]*vector_2[2]);
  double dot = PDM_DOT_PRODUCT(vector_1,vector_2);
  double cross[3];

  if (PDM_ABS(dot) > ((norm_1*norm_2)*(1.-QUATERNION_EPS))) {
    if (dot > 0.) { // parallel vectors
      PDM_quaternion_set_identity(qt_out);
      return;
    }
    else { // opposite vectors
      double z_axis[3] = {0., 0., 1.};
      double dot2 = PDM_DOT_PRODUCT(vector_1,z_axis);
      if (PDM_ABS(dot2) > (norm_1*(1-QUATERNION_EPS))) { // vector 1 is z_axis
        double y_axis[3] = {0., 1., 0.};
        PDM_CROSS_PRODUCT(cross,vector_1,y_axis);
      }
      else { // vector != z_axis -> rotation around z-axis
        PDM_CROSS_PRODUCT(cross,vector_1,z_axis);
      }
      qt_out->q[0] = 0.;
      qt_out->q_squared[0] = 0.;
    }
  }
  else {
    PDM_CROSS_PRODUCT(cross,vector_1,vector_2);
    qt_out->q[0] = sqrt(normsq_1*normsq_2)+dot;
    qt_out->q_squared[0] = normsq_1*normsq_2+dot*dot+2*sqrt(normsq_1*normsq_2)*dot;
  }
  qt_out->q[1] = cross[0];
  qt_out->q[2] = cross[1];
  qt_out->q[3] = cross[2];
  qt_out->q_squared[1] = cross[0]*cross[0];
  qt_out->q_squared[2] = cross[1]*cross[1];
  qt_out->q_squared[3] = cross[2]*cross[2];
  PDM_quaternion_normalize(qt_out);
}

void
PDM_quaternion_from_axis_angle
(
  const double          axis[3],
  const double          angle,
        PDM_quaternion* qt_out
)
{
  qt_out->q[0] = cos(0.5*angle);
  qt_out->q_squared[0] = qt_out->q[0] * qt_out->q[0];
  double c         = sin(0.5*angle);
  double c_squared = c*c;
  double axis_invnorm = 1./sqrt(PDM_DOT_PRODUCT(axis,axis));
  double axis_invnorm_sq = 1./PDM_DOT_PRODUCT(axis,axis);
  qt_out->q[1] = c*axis[0]*axis_invnorm;
  qt_out->q[2] = c*axis[1]*axis_invnorm;
  qt_out->q[3] = c*axis[2]*axis_invnorm;
  qt_out->q_squared[1] = c_squared*axis[0]*axis[0]*axis_invnorm_sq;
  qt_out->q_squared[2] = c_squared*axis[1]*axis[1]*axis_invnorm_sq;
  qt_out->q_squared[3] = c_squared*axis[2]*axis[2]*axis_invnorm_sq;
}

void
PDM_quaternion_from_axis_angle_derivative(
  const double          axis[3],
  const double          angle,
  const double          axis_der[3],
  const double          angle_der,
        PDM_quaternion* qt_out
)
{
  // /!\ not tested in pdm_quaternion.test.cpp yet !
  double c     =      sin(0.5*angle);
  double c_der =  0.5*cos(0.5*angle)*angle_der;
  qt_out->q[0]         = -0.5*sin(0.5*angle)*angle_der;

  double ax_sq_norm = PDM_DOT_PRODUCT(axis,axis);
  double axis_invnorm     = 1./sqrt(ax_sq_norm);
  double axis_invnorm_der = -(axis[0]*axis_der[0]+axis[1]*axis_der[1]+axis[2]*axis_der[2])/(ax_sq_norm*sqrt(ax_sq_norm));
  qt_out->q[1] = c_der*axis[0]*axis_invnorm + c*axis_der[0]*axis_invnorm+c*axis[0]*axis_invnorm_der;
  qt_out->q[2] = c_der*axis[1]*axis_invnorm + c*axis_der[1]*axis_invnorm+c*axis[1]*axis_invnorm_der;
  qt_out->q[3] = c_der*axis[2]*axis_invnorm + c*axis_der[2]*axis_invnorm+c*axis[2]*axis_invnorm_der;

  // TODO: improve precision of squared part
  _PDM_quaternion_compute_squared(qt_out);
}

void
PDM_quaternion_from_axis_aligned_rotation
(
  const double angle,
  const int axis_ind,
  PDM_quaternion* qt_out
)
{
  qt_out->q[0] = cos(angle*0.5);
  qt_out->q[1] = 0.;
  qt_out->q[2] = 0.;
  qt_out->q[3] = 0.;
  qt_out->q[axis_ind+1] = sin(angle*0.5);

  _PDM_quaternion_compute_squared(qt_out);
}

void
PDM_quaternion_from_x_rotation
(
  const double          angle,
        PDM_quaternion* qt
)
{
  PDM_quaternion_from_axis_aligned_rotation(angle,0,qt);
}

void
PDM_quaternion_from_y_rotation
(
  const double          angle,
        PDM_quaternion* qt
)
{
  PDM_quaternion_from_axis_aligned_rotation(angle,1,qt);
}

void
PDM_quaternion_from_z_rotation
(
  const double          angle,
        PDM_quaternion* qt
)
{
  PDM_quaternion_from_axis_aligned_rotation(angle,2,qt);
}

void
PDM_quaternion_from_axis_aligned_symmetry
(
  const int             axis_ind,
        PDM_quaternion* qt
)
{
  qt->q[0] = 0.;
  qt->q[1] = 0.;
  qt->q[2] = 0.;
  qt->q[3] = 0.;
  qt->q[axis_ind+1] = 1.;

  _PDM_quaternion_compute_squared(qt);
}

void
PDM_quaternion_from_x_symmetry
(
  PDM_quaternion* qt
)
{
  PDM_quaternion_from_axis_aligned_symmetry(0,qt);
}

void
PDM_quaternion_from_y_symmetry
(
  PDM_quaternion* qt
)
{
  PDM_quaternion_from_axis_aligned_symmetry(1,qt);
}

void
PDM_quaternion_from_z_symmetry
(
  PDM_quaternion* qt
)
{
  PDM_quaternion_from_axis_aligned_symmetry(2,qt);
}

void
PDM_quaternion_from_euler_angles
(
  const double          ang_x,
  const double          ang_y,
  const double          ang_z,
  const int             order[3],
  const PDM_bool_t      intrinsic,
        PDM_quaternion* qt
)
{
  double angles[3] = {ang_x,ang_y,ang_z};
  PDM_quaternion_set_identity(qt);
  PDM_quaternion q_curr;
  for (size_t i = 0; i < 3; i++) {
    PDM_quaternion_from_axis_aligned_rotation(angles[order[i]],order[i],&q_curr);
    if ( intrinsic ) {
      PDM_quaternion_compose(qt,&q_curr,qt);
    }else {
      PDM_quaternion_compose(&q_curr,qt,qt);
    }
  }
}

void
PDM_quaternion_from_rotation_matrix
(
  const double*         rotation_matrix,
        PDM_quaternion* qt
)
{
  double tr = rotation_matrix[3*0+0]+rotation_matrix[3*1+1]+rotation_matrix[3*2+2];
  // printf("%s::%d det = %23.16e\n",__FILE__,__LINE__,_PDM_quaternion_3x3_mat_det(rotation_matrix));
  double pivots[4] = {
    rotation_matrix[3*0+0],
    rotation_matrix[3*1+1],
    rotation_matrix[3*2+2],
    tr,
  };
  int choice = 0;
  for (int i=1;i<4;i++){
    if (pivots[i]>pivots[choice]){
      choice = i;
    }
  };

  int i,j,k;

  if (choice != 3){ // 0,1,2
    i = choice;
    j = (i+1)%3;
    k = (j+1)%3;

    qt->q[0]   = rotation_matrix[3*k+j] - rotation_matrix[3*j+k]; // w
    qt->q[1+i] = 1 + rotation_matrix[3*i+i] - rotation_matrix[3*j+j] - rotation_matrix[3*k+k];
    qt->q[1+j] = rotation_matrix[3*j+i] + rotation_matrix[3*i+j];
    qt->q[1+k] = rotation_matrix[3*k+i] + rotation_matrix[3*i+k];
  }
  else {
    qt->q[0] = 1+tr; // w
    qt->q[1] = rotation_matrix[3*2+1] - rotation_matrix[3*1+2];
    qt->q[2] = rotation_matrix[3*0+2] - rotation_matrix[3*2+0];
    qt->q[3] = rotation_matrix[3*1+0] - rotation_matrix[3*0+1];
  }
  _PDM_quaternion_compute_squared(qt);
  PDM_quaternion_normalize(qt);

}

void
PDM_quaternion_from_homogeneous_matrix
(
  const double*         homogeneous_matrix,
        PDM_quaternion* qt
)
{
  double rot_matrix[9];
  rot_matrix[3*0+0] = homogeneous_matrix[4*0+0];
  rot_matrix[3*0+1] = homogeneous_matrix[4*0+1];
  rot_matrix[3*0+2] = homogeneous_matrix[4*0+2];
  rot_matrix[3*1+0] = homogeneous_matrix[4*1+0];
  rot_matrix[3*1+1] = homogeneous_matrix[4*1+1];
  rot_matrix[3*1+2] = homogeneous_matrix[4*1+2];
  rot_matrix[3*2+0] = homogeneous_matrix[4*2+0];
  rot_matrix[3*2+1] = homogeneous_matrix[4*2+1];
  rot_matrix[3*2+2] = homogeneous_matrix[4*2+2];
  PDM_quaternion_from_rotation_matrix(rot_matrix,qt);

}

/*----------------------------------------------------------------------------
 *  PDM_quaternion_t TO OTHER FORMATS
 *----------------------------------------------------------------------------*/

void
PDM_quaternion_to_axis_angle
(
  const PDM_quaternion* qt,
        double  axis[3],
        double* angle
)
{
  // Prefer positive angle
  (*angle) = 2.0 * acos(qt->q[0]);
  double sign = (*angle)>0. ? 1. : -1.;
  (*angle) *= sign;
  double divider = sign*sqrt(1.0 - qt->q_squared[0]);

  if(PDM_ABS(divider) > QUATERNION_EPS) {
    // Calculate the axis
    axis[0] = qt->q[1] / divider;
    axis[1] = qt->q[2] / divider;
    axis[2] = qt->q[3] / divider;
  } else {
    // Arbitrary normalized axis
    axis[0] = 1.;
    axis[1] = 0.;
    axis[2] = 0.;
  }
}

void
PDM_quaternion_to_euler_angles
(
  const PDM_quaternion* qt,
  const int             order[3],
  const PDM_bool_t      intrinsic,
        double*         ang_x,
        double*         ang_y,
        double*         ang_z
)
{
  int i = intrinsic ? order[2] : order[0];
  int j = order[1];
  int k = intrinsic ? order[0] : order[2];
  PDM_bool_t symmetric = (PDM_bool_t)(i == k);
  if ( symmetric ) {
      k = 3-i-j;
  }
  double sign = (i-j)*(j-k)*(k-i)/2;
  double quat[4] = {qt->q[1],qt->q[2],qt->q[3],qt->q[0]};
  double a,b,c,d;
  if (symmetric){
      a = quat[3];
      b = quat[i];
      c = quat[j];
      d = quat[k];
  }else{
      a = quat[3]-quat[j];
      b = quat[i]+quat[k]*sign;
      c = quat[j]+quat[3];
      d = quat[k]*sign-quat[i];
  }
  double angle[3];
  double half_sum = atan2(b,a);
  double half_diff = atan2(d,c);

  // second angle
  angle[1] = 2*atan2(hypot(c,d),hypot(a,b));
  if (PDM_ABS(angle[1])<=1e-7){
      angle[2] = 0.;
      angle[0] = 2*half_sum;
  }else if (PDM_ABS(angle[1]-PDM_PI)<=1e-7){
      angle[2] = 0.;
      angle[0] = 2*half_diff*(intrinsic?1.:-1.);
      /* code */
  }else{
      if (intrinsic){
          angle[2] = half_sum-half_diff;
          angle[0] = half_sum+half_diff;
      }else{
          angle[0] = half_sum-half_diff;
          angle[2] = half_sum+half_diff;

      }
  }
  if (!symmetric){
      angle[1] -= 0.5*PDM_PI;
      if (intrinsic){
          angle[0] *= sign;
      }else{
          angle[2] *= sign;
      }
  }
  if (angle[0] < -PDM_PI) angle[0] += 2*PDM_PI;
  if (angle[0] >  PDM_PI) angle[0] -= 2*PDM_PI;
  if (angle[1] < -PDM_PI) angle[1] += 2*PDM_PI;
  if (angle[1] >  PDM_PI) angle[1] -= 2*PDM_PI;
  if (angle[2] < -PDM_PI) angle[2] += 2*PDM_PI;
  if (angle[2] >  PDM_PI) angle[2] -= 2*PDM_PI;

  int rev_order[3];
  rev_order[order[0]] = 0;
  rev_order[order[1]] = 1;
  rev_order[order[2]] = 2;

  (*ang_x) = angle[rev_order[0]];
  (*ang_y) = angle[rev_order[1]];
  (*ang_z) = angle[rev_order[2]];

}

void
PDM_quaternion_to_rotation_matrix
(
  const PDM_quaternion* qt,
        double*         rotation_matrix
)
{
  rotation_matrix[3*0+0] = 1 - 2 * qt->q_squared[2]    - 2 * qt->q_squared[3]   ;
  rotation_matrix[3*0+1] =     2 * qt->q[1] * qt->q[2] - 2 * qt->q[3] * qt->q[0];
  rotation_matrix[3*0+2] =     2 * qt->q[1] * qt->q[3] + 2 * qt->q[2] * qt->q[0];
  rotation_matrix[3*1+0] =     2 * qt->q[1] * qt->q[2] + 2 * qt->q[3] * qt->q[0];
  rotation_matrix[3*1+1] = 1 - 2 * qt->q_squared[3]    - 2 * qt->q_squared[1]   ;
  rotation_matrix[3*1+2] =     2 * qt->q[2] * qt->q[3] - 2 * qt->q[1] * qt->q[0];
  rotation_matrix[3*2+0] =     2 * qt->q[1] * qt->q[3] - 2 * qt->q[2] * qt->q[0];
  rotation_matrix[3*2+1] =     2 * qt->q[2] * qt->q[3] + 2 * qt->q[1] * qt->q[0];
  rotation_matrix[3*2+2] = 1 - 2 * qt->q_squared[1]    - 2 * qt->q_squared[2]   ;
}

void
PDM_quaternion_to_homogeneous_matrix
(
  const PDM_quaternion* qt,
        double*         homogeneous_matrix
)
{
  double rot_mat[9];
  PDM_quaternion_to_rotation_matrix(qt,rot_mat);
  homogeneous_matrix[4*0+0] = rot_mat[3*0+0];
  homogeneous_matrix[4*0+1] = rot_mat[3*0+1];
  homogeneous_matrix[4*0+2] = rot_mat[3*0+2];
  homogeneous_matrix[4*0+3] = 0.;
  homogeneous_matrix[4*1+0] = rot_mat[3*1+0];
  homogeneous_matrix[4*1+1] = rot_mat[3*1+1];
  homogeneous_matrix[4*1+2] = rot_mat[3*1+2];
  homogeneous_matrix[4*1+3] = 0.;
  homogeneous_matrix[4*2+0] = rot_mat[3*2+0];
  homogeneous_matrix[4*2+1] = rot_mat[3*2+1];
  homogeneous_matrix[4*2+2] = rot_mat[3*2+2];
  homogeneous_matrix[4*2+3] = 0.;
  homogeneous_matrix[4*3+0] = 0.;
  homogeneous_matrix[4*3+1] = 0.;
  homogeneous_matrix[4*3+2] = 0.;
  homogeneous_matrix[4*3+3] = 1.;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
