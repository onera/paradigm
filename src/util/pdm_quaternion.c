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
        double qt[4]
)
{
  qt[0] = w;
  qt[1] = v1;
  qt[2] = v2;
  qt[3] = v3;
}

void
PDM_quaternion_set_identity
(
  double qt[4]
)
{
  qt[0] = 1.;
  qt[1] = 0.;
  qt[2] = 0.;
  qt[3] = 0.;
}

void
PDM_quaternion_print
(
  const double qt[4]
)
{
  PDM_printf("qt = [%12.5e %12.5e %12.5e %12.5e]\n",
    qt[0],
    qt[1],
    qt[2],
    qt[3]);
}

void
PDM_quaternion_conjugate
(
  const double qt[4],
        double qt_out[4]
)
{
  qt_out[0] =  qt[0];
  qt_out[1] = -qt[1];
  qt_out[2] = -qt[2];
  qt_out[3] = -qt[3];
}

double
PDM_quaternion_norm
(
  const double qt[4]
)
{
  return sqrt(
    qt[0] * qt[0]  +
    qt[1] * qt[1] + 
    qt[2] * qt[2] + 
    qt[3] * qt[3]);
}

void
PDM_quaternion_normalize
(
  double qt[4]
)
{
  double len = PDM_quaternion_norm(qt);
  qt[0] = qt[0] / len;
  qt[1] = qt[1] / len;
  qt[2] = qt[2] / len;
  qt[3] = qt[3] / len;
}

void
PDM_quaternion_compose
(
  const double qt_1[4],
  const double qt_2[4],
        double qt_out[4]
)
{
  // buffering data in case where qt_1==qt_out or qt_2==qt_out
  double w,v1,v2,v3;
  w  = qt_1[0] * qt_2[0] - qt_1[1] * qt_2[1] - qt_1[2] * qt_2[2] - qt_1[3] * qt_2[3];
  v1 = qt_1[1] * qt_2[0] + qt_1[0] * qt_2[1] + qt_1[2] * qt_2[3] - qt_1[3] * qt_2[2];
  v2 = qt_1[0] * qt_2[2] - qt_1[1] * qt_2[3] + qt_1[2] * qt_2[0] + qt_1[3] * qt_2[1];
  v3 = qt_1[0] * qt_2[3] + qt_1[1] * qt_2[2] - qt_1[2] * qt_2[1] + qt_1[3] * qt_2[0] ;
  qt_out[0] = w ;
  qt_out[1] = v1;
  qt_out[2] = v2;
  qt_out[3] = v3;
}

PDM_bool_t
PDM_quaternion_equal
(
  const double qt[4],
  const double w,
  const double v0,
  const double v1,
  const double v2,
  const double epsilon
)
{
  return (PDM_ABS(qt[0] - w ) < epsilon) & \
         (PDM_ABS(qt[1] - v0) < epsilon) &\
         (PDM_ABS(qt[2] - v1) < epsilon) &\
         (PDM_ABS(qt[3] - v2) < epsilon);
}

PDM_bool_t
PDM_quaternion_equal_quaternion
(
  const double qt_1[4],
  const double qt_2[4],
  const double epsilon
)
{
  return PDM_quaternion_equal(qt_1,qt_2[0],qt_2[1],qt_2[2],qt_2[3],epsilon);
}


/*----------------------------------------------------------------------------
 *  APPLICATION TO VECTORS
 *----------------------------------------------------------------------------*/

void
PDM_quaternion_rotate
(
  const double qt[4],
  const double* vector,
  const int n_samp,
        double* vector_out
)
{
  double ww = qt[0] * qt[0];
  double xx = qt[1] * qt[1];
  double yy = qt[2] * qt[2];
  double zz = qt[3] * qt[3];
  double wx = qt[0] * qt[1];
  double wy = qt[0] * qt[2];
  double wz = qt[0] * qt[3];
  double xy = qt[1] * qt[2];
  double xz = qt[1] * qt[3];
  double yz = qt[2] * qt[3];

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
  const double qt[4],
  const double qt_der[4],
  const double* vector,
  const int n_samp,
        double* vector_output_der
)
{
  double ww_der = 2*qt[0] * qt_der[0];
  double xx_der = 2*qt[1] * qt_der[1];
  double yy_der = 2*qt[2] * qt_der[2];
  double zz_der = 2*qt[3] * qt_der[3];
  double wx_der = qt_der[0] * qt[1] + qt[0] * qt_der[1];
  double wy_der = qt_der[0] * qt[2] + qt[0] * qt_der[2];
  double wz_der = qt_der[0] * qt[3] + qt[0] * qt_der[3];
  double xy_der = qt_der[1] * qt[2] + qt[1] * qt_der[2];
  double xz_der = qt_der[1] * qt[3] + qt[1] * qt_der[3];
  double yz_der = qt_der[2] * qt[3] + qt[2] * qt_der[3];

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
  double qt[4];
  PDM_quaternion_from_two_vectors(vector_1,vector_2,qt);
  double axis[3];
  double angle;
  PDM_quaternion_to_axis_angle(qt,axis,&angle);
  for (int i = 0; i < n_samp; i++) {
    PDM_quaternion_from_axis_angle(axis,t[i]*angle,qt);
    PDM_quaternion_rotate(qt,vector_1,1,&out[3*i]);
  }
}

void
PDM_quaternion_slerp_from_two_vectors_derivative
(
  const double vector_1[3],
  const double vector_2[3],
  const double* t,
  const int n_samp,
        double* out
)
{
  double qt[4];
  double qt_der[4];
  PDM_quaternion_from_two_vectors(vector_1,vector_2,qt);
  double axis[3];
  double axis_der[3] = {0,0,0};
  double angle;
  PDM_quaternion_to_axis_angle(qt,axis,&angle);
  for (int i = 0; i < n_samp; i++) {
    PDM_quaternion_from_axis_angle(axis,t[i]*angle,qt);
    PDM_quaternion_from_axis_angle_derivative(axis,t[i]*angle,axis_der,angle,qt_der);
    PDM_quaternion_rotate_derivative(qt,qt_der,vector_1,1,&out[3*i]);
  } 
}


/*----------------------------------------------------------------------------
 *  PDM_quaternion_t FROM OTHER FORMATS
 *----------------------------------------------------------------------------*/

void 
PDM_quaternion_from_two_vectors
(
  const double* vector_1,
  const double* vector_2,
        double  qt_out[4]
)
{
  double norm_1 = sqrt(vector_1[0]*vector_1[0]+vector_1[1]*vector_1[1]+vector_1[2]*vector_1[2]);
  double norm_2 = sqrt(vector_2[0]*vector_2[0]+vector_2[1]*vector_2[1]+vector_2[2]*vector_2[2]);
  double norm_vector_1[3];
  double norm_vector_2[3];
  for (int i = 0; i < 3; i++) {
    norm_vector_1[i] = vector_1[i]/norm_1;
    norm_vector_2[i] = vector_2[i]/norm_2;
  }
  double dot = PDM_DOT_PRODUCT(norm_vector_1,norm_vector_2);
  double cross[3];
  PDM_CROSS_PRODUCT(cross,norm_vector_1,norm_vector_2);

  if (PDM_ABS(dot) > 1.-QUATERNION_EPS) {
    if (dot > 0.) { // parallel vectors
      PDM_quaternion_set_identity(qt_out);
    }
    else { // opposite vectors
      double y_axis[3] = {0., 1., 0.};
      double dot2 = PDM_DOT_PRODUCT(norm_vector_1,y_axis);
      if (PDM_ABS(dot2) > 1. - QUATERNION_EPS) { // vector 1 is y_axis
          double axis[3] = {1.,0.,0.}; // x_axis
          PDM_quaternion_from_axis_angle(axis,M_PI,qt_out);
      }
      else { // vector != y_axis
        double axis[3];
        PDM_CROSS_PRODUCT(axis,y_axis,norm_vector_1);
        PDM_quaternion_from_axis_angle(axis,M_PI,qt_out);
      }
    }
  }
  else {
    qt_out[0] = 1+dot;
    qt_out[1] = cross[0];
    qt_out[2] = cross[1];
    qt_out[3] = cross[2];
    PDM_quaternion_normalize(qt_out);
  }
}

void 
PDM_quaternion_from_axis_angle
(
  const double axis[3],
  const double angle,
        double qt_out[4]
)
{ 
  qt_out[0] = cos(0.5*angle);
  double c  = sin(0.5*angle);
  double axis_invnorm = 1./sqrt(PDM_DOT_PRODUCT(axis,axis));
  qt_out[1] = c * axis[0]*axis_invnorm;
  qt_out[2] = c * axis[1]*axis_invnorm;
  qt_out[3] = c * axis[2]*axis_invnorm;
}

void
PDM_quaternion_from_axis_angle_derivative(
  const double axis[3], 
  const double angle,
  const double axis_der[3],
  const double angle_der,
        double qt_out[4]
)
{
  double c     =      sin(0.5*angle);
  qt_out[0]    = -0.5*sin(0.5*angle)*angle_der;
  double c_der =  0.5*cos(0.5*angle)*angle_der;
  double ax_sq_norm = PDM_DOT_PRODUCT(axis,axis);
  double axis_invnorm     = 1./sqrt(ax_sq_norm);
  double axis_invnorm_der = -(axis[0]*axis_der[0]+axis[1]*axis_der[1]+axis[2]*axis_der[2])/(ax_sq_norm*sqrt(ax_sq_norm));
  qt_out[1] = c_der*axis[0]*axis_invnorm + c*axis_der[0]*axis_invnorm+c*axis[0]*axis_invnorm_der;
  qt_out[2] = c_der*axis[1]*axis_invnorm + c*axis_der[1]*axis_invnorm+c*axis[1]*axis_invnorm_der;
  qt_out[3] = c_der*axis[2]*axis_invnorm + c*axis_der[2]*axis_invnorm+c*axis[2]*axis_invnorm_der;
}

void
PDM_quaternion_from_axis_aligned_rotation
(
  const double angle,
  const int axis_ind,
        double qt_out[4]
)
{
  qt_out[0]  = cos(angle*0.5);
  qt_out[1] = 0.;
  qt_out[2] = 0.;
  qt_out[3] = 0.;
  qt_out[axis_ind+1] = sin(angle*0.5);
}

void
PDM_quaternion_from_x_rotation
(
  const double angle,
        double qt[4]
)
{
  PDM_quaternion_from_axis_aligned_rotation(angle,0,qt);
}

void
PDM_quaternion_from_y_rotation
(
  const double angle,
        double qt[4]
)
{
  PDM_quaternion_from_axis_aligned_rotation(angle,1,qt);
}

void
PDM_quaternion_from_z_rotation
(
  const double angle,
        double qt[4]
)
{
  PDM_quaternion_from_axis_aligned_rotation(angle,2,qt);
}

void
PDM_quaternion_from_axis_aligned_symmetry
(
  const int axis_ind,
        double qt[4]
)
{
  qt[0] = 0.;
  qt[1] = 0.;
  qt[2] = 0.;
  qt[3] = 0.;
  qt[axis_ind+1] = 1.;
}

void
PDM_quaternion_from_x_symmetry
(
  double qt[4]
)
{
  PDM_quaternion_from_axis_aligned_symmetry(0,qt);
}

void
PDM_quaternion_from_y_symmetry
(
  double qt[4]
)
{
  PDM_quaternion_from_axis_aligned_symmetry(1,qt);
}

void
PDM_quaternion_from_z_symmetry
(
  double qt[4]
)
{
  PDM_quaternion_from_axis_aligned_symmetry(2,qt);
}

void 
PDM_quaternion_from_euler_angles
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  PDM_bool_t intrinsic,
        double qt[4]
)
{
  double angles[3] = {ang_x,ang_y,ang_z};
  PDM_quaternion_set_identity(qt);
  double q_curr[4];
  for (size_t i = 0; i < 3; i++) {
    PDM_quaternion_from_axis_aligned_rotation(angles[order[i]],order[i],q_curr);
    if ( intrinsic ) {
      PDM_quaternion_compose(qt,q_curr,qt);
    }else {
      PDM_quaternion_compose(q_curr,qt,qt);
    }
  }
}

void
PDM_quaternion_from_rotation_matrix
(
  const double* rotation_matrix,
        double qt_out[4]
)
{
  double tr = rotation_matrix[3*0+0]+rotation_matrix[3*1+1]+rotation_matrix[3*2+2];
  double tmp;
  if (PDM_ABS(tr) > QUATERNION_EPS){
    tmp = sqrt(1+tr)*2.;
    qt_out[0]    = .25*tmp;
    qt_out[1] = (rotation_matrix[3*2+1]-rotation_matrix[3*1+2])/tmp;
    qt_out[2] = (rotation_matrix[3*0+2]-rotation_matrix[3*2+0])/tmp;
    qt_out[3] = (rotation_matrix[3*1+0]-rotation_matrix[3*0+1])/tmp;
  }
  else if ((rotation_matrix[3*0+0]>rotation_matrix[3*1+1]) & (rotation_matrix[3*0+0]>rotation_matrix[3*2+2])){
    tmp = sqrt(1+rotation_matrix[3*0+0]-rotation_matrix[3*1+1]-rotation_matrix[3*2+2])*2.;
    qt_out[0]    = (rotation_matrix[3*2+1]-rotation_matrix[3*1+2])/tmp;
    qt_out[1] = .25*tmp;
    qt_out[2] = (rotation_matrix[3*1+0]+rotation_matrix[3*0+1])/tmp;
    qt_out[3] = (rotation_matrix[3*2+0]+rotation_matrix[3*0+2])/tmp;
  }
  else if (rotation_matrix[3*1+1]>rotation_matrix[3*2+2]){
    tmp  = sqrt(1+rotation_matrix[3*1+1] - rotation_matrix[3*0+0] - rotation_matrix[3*2+2]) * 2;
    qt_out[0]    = (rotation_matrix[3*0+2] - rotation_matrix[3*2+0]) / tmp;
    qt_out[1] = (rotation_matrix[3*0+1] + rotation_matrix[3*1+0]) / tmp; 
    qt_out[2] = 0.25 * tmp;
    qt_out[3] = (rotation_matrix[3*1+2] + rotation_matrix[3*2+1]) / tmp; 
  } else { 
    tmp = sqrt(1.0 + rotation_matrix[3*2+2] - rotation_matrix[3*0+0] - rotation_matrix[3*1+1]) * 2;
    qt_out[0]    = (rotation_matrix[3*1+0] - rotation_matrix[3*0+1]) / tmp;
    qt_out[1] = (rotation_matrix[3*0+2] + rotation_matrix[3*2+0]) / tmp;
    qt_out[2] = (rotation_matrix[3*1+2] + rotation_matrix[3*2+1]) / tmp;
    qt_out[3] = 0.25 * tmp;
  }
}

void
PDM_quaternion_from_homogeneous_matrix
(
  const double *homogeneous_matrix,
        double qt_out[4]
)
{
  double tr = homogeneous_matrix[4*0+0]+homogeneous_matrix[4*1+1]+homogeneous_matrix[4*2+2];
  double tmp;
  if (PDM_ABS(tr) > QUATERNION_EPS){
    tmp = sqrt(1+tr)*2.;
    qt_out[0]    = .25*tmp;
    qt_out[1] = (homogeneous_matrix[4*2+1]-homogeneous_matrix[4*1+2])/tmp;
    qt_out[2] = (homogeneous_matrix[4*0+2]-homogeneous_matrix[4*2+0])/tmp;
    qt_out[3] = (homogeneous_matrix[4*1+0]-homogeneous_matrix[4*0+1])/tmp;
  }
  else if ((homogeneous_matrix[4*0+0]>homogeneous_matrix[4*1+1]) & (homogeneous_matrix[4*0+0]>homogeneous_matrix[4*2+2])){
    tmp = sqrt(1+homogeneous_matrix[4*0+0]-homogeneous_matrix[4*1+1]-homogeneous_matrix[4*2+2])*2.;
    qt_out[0]    = (homogeneous_matrix[4*2+1]-homogeneous_matrix[4*1+2])/tmp;
    qt_out[1] = .25*tmp;
    qt_out[2] = (homogeneous_matrix[4*1+0]+homogeneous_matrix[4*0+1])/tmp;
    qt_out[3] = (homogeneous_matrix[4*2+0]+homogeneous_matrix[4*0+2])/tmp;
  }
  else if (homogeneous_matrix[4*1+1]>homogeneous_matrix[4*2+2]){
    tmp  = sqrt(1+homogeneous_matrix[4*1+1] - homogeneous_matrix[4*0+0] - homogeneous_matrix[4*2+2]) * 2;
    qt_out[0]    = (homogeneous_matrix[4*0+2] - homogeneous_matrix[4*2+0]) / tmp;
    qt_out[1] = (homogeneous_matrix[4*0+1] + homogeneous_matrix[4*1+0]) / tmp; 
    qt_out[2] = 0.25 * tmp;
    qt_out[3] = (homogeneous_matrix[4*1+2] + homogeneous_matrix[4*2+1]) / tmp; 
  } else { 
    tmp = sqrt(1.0 + homogeneous_matrix[4*2+2] - homogeneous_matrix[4*0+0] - homogeneous_matrix[4*1+1]) * 2;
    qt_out[0]    = (homogeneous_matrix[4*1+0] - homogeneous_matrix[4*0+1]) / tmp;
    qt_out[1] = (homogeneous_matrix[4*0+2] + homogeneous_matrix[4*2+0]) / tmp;
    qt_out[2] = (homogeneous_matrix[4*1+2] + homogeneous_matrix[4*2+1]) / tmp;
    qt_out[3] = 0.25 * tmp;
  }
}

/*----------------------------------------------------------------------------
 *  PDM_quaternion_t TO OTHER FORMATS
 *----------------------------------------------------------------------------*/

void
PDM_quaternion_to_axis_angle
(
  const double qt[4],
        double  axis[3],
        double* angle
)
{
  (*angle) = 2.0 * acos(qt[0]);
  double divider = sqrt(1.0 - qt[0] * qt[0]);

  if(PDM_ABS(divider) > QUATERNION_EPS) {
    // Calculate the axis
    axis[0] = qt[1] / divider;
    axis[1] = qt[2] / divider;
    axis[2] = qt[3] / divider;
  } else {
    // Arbitrary normalized axis
    axis[0] = 1;
    axis[1] = 0;
    axis[2] = 0;
  }
}

void
PDM_quaternion_to_euler_angles
(
  const double qt[4],
  const int order[3],
  const PDM_bool_t intrinsic,
        double* ang_x,
        double* ang_y,
        double* ang_z
)
{
    int i = intrinsic ? order[2] : order[0];
    int j = order[1];
    int k = intrinsic ? order[0] : order[2];
    PDM_bool_t symmetric = i==k;
    if ( symmetric ) {
        k = 3-i-j;
    }
    double sign = (i-j)*(j-k)*(k-i)/2;
    double quat[4] = {qt[1],qt[2],qt[3],qt[0]};
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
    }else if (PDM_ABS(angle[1]-M_PI)<=1e-7){
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
        angle[1] -= 0.5*M_PI;
        if (intrinsic){
            angle[0] *= sign;
        }else{
            angle[2] *= sign;
        }
    }
    if (angle[0] < -M_PI) angle[0] += 2*M_PI;
    if (angle[0] >  M_PI) angle[0] -= 2*M_PI;
    if (angle[1] < -M_PI) angle[1] += 2*M_PI;
    if (angle[1] >  M_PI) angle[1] -= 2*M_PI;
    if (angle[2] < -M_PI) angle[2] += 2*M_PI;
    if (angle[2] >  M_PI) angle[2] -= 2*M_PI;

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
  const double qt[4], 
        double* rotation_matrix
)
{
  rotation_matrix[3*0+0] = 1 - 2 * qt[2] * qt[2] - 2 * qt[3] * qt[3];
  rotation_matrix[3*0+1] =     2 * qt[1] * qt[2] - 2 * qt[3] * qt[0] ;
  rotation_matrix[3*0+2] =     2 * qt[1] * qt[3] + 2 * qt[2] * qt[0] ;
  rotation_matrix[3*1+0] =     2 * qt[1] * qt[2] + 2 * qt[3] * qt[0] ;
  rotation_matrix[3*1+1] = 1 - 2 * qt[3] * qt[3] - 2 * qt[1] * qt[1];
  rotation_matrix[3*1+2] =     2 * qt[2] * qt[3] - 2 * qt[1] * qt[0] ;
  rotation_matrix[3*2+0] =     2 * qt[1] * qt[3] - 2 * qt[2] * qt[0] ;
  rotation_matrix[3*2+1] =     2 * qt[2] * qt[3] + 2 * qt[1] * qt[0] ;
  rotation_matrix[3*2+2] = 1 - 2 * qt[1] * qt[1] - 2 * qt[2] * qt[2];
}

void
PDM_quaternion_to_homogeneous_matrix
(
  const double qt[4],
        double* homogeneous_matrix
)
{
  homogeneous_matrix[4*0+0] = 1 - 2 * qt[2] * qt[2] - 2 * qt[3] * qt[3];
  homogeneous_matrix[4*0+1] =     2 * qt[1] * qt[2] - 2 * qt[3] * qt[0] ;
  homogeneous_matrix[4*0+2] =     2 * qt[1] * qt[3] + 2 * qt[2] * qt[0] ;
  homogeneous_matrix[4*0+3] = 0.;
  homogeneous_matrix[4*1+0] =     2 * qt[1] * qt[2] + 2 * qt[3] * qt[0] ;
  homogeneous_matrix[4*1+1] = 1 - 2 * qt[3] * qt[3] - 2 * qt[1] * qt[1];
  homogeneous_matrix[4*1+2] =     2 * qt[2] * qt[3] - 2 * qt[1] * qt[0] ;
  homogeneous_matrix[4*1+3] = 0.;
  homogeneous_matrix[4*2+0] =     2 * qt[1] * qt[3] - 2 * qt[2] * qt[0] ;
  homogeneous_matrix[4*2+1] =     2 * qt[2] * qt[3] + 2 * qt[1] * qt[0] ;
  homogeneous_matrix[4*2+2] = 1 - 2 * qt[1] * qt[1] - 2 * qt[2] * qt[2];
  homogeneous_matrix[4*2+3] = 0.;
  homogeneous_matrix[4*3+0] = 0.;
  homogeneous_matrix[4*3+1] = 0.;
  homogeneous_matrix[4*3+2] = 0.;
  homogeneous_matrix[4*3+3] = 1.;
}



#ifdef __cplusplus
}
#endif /* __cplusplus */
