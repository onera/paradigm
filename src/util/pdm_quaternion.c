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




/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_quaternion_t*
PDM_quaternion_create
(
  const double w,
  const double v1,
  const double v2,
  const double v3
)
{
  PDM_quaternion_t *qt;
  PDM_malloc(qt, 1, PDM_quaternion_t);

  qt->w  = w;
  qt->v[0] = v1;
  qt->v[1] = v2;
  qt->v[2] = v3;

  return qt;
}

void
PDM_quaternion_free
(
  PDM_quaternion_t* qt
)
{
  if (qt == NULL) {
    return;
  }

  free(qt);
}

void
PDM_quaternion_set_identity
(
  PDM_quaternion_t* qt
)
{
  qt->w  = 1.;
  qt->v[0] = 0.;
  qt->v[1] = 0.;
  qt->v[2] = 0.;
}

void
PDM_quaternion_print
(
  const PDM_quaternion_t* qt
)
{
  PDM_printf("qt = [%12.5e %12.5e %12.5e %12.5e]\n",
    qt->w   ,
    qt->v[0],
    qt->v[1],
    qt->v[2]);
}

void
PDM_quaternion_conjugate
(
  const PDM_quaternion_t* qt,
  PDM_quaternion_t* qt_out
)
{
  qt_out->w  =  qt->w;
  qt_out->v[0] = -qt->v[0];
  qt_out->v[1] = -qt->v[1];
  qt_out->v[2] = -qt->v[2];
}

double
PDM_quaternion_norm
(
  const PDM_quaternion_t* qt
)
{
  return sqrt(
    qt->w  * qt->w  +
    qt->v[0] * qt->v[0] + 
    qt->v[1] * qt->v[1] + 
    qt->v[2] * qt->v[2]);
}

void
PDM_quaternion_normalize
(
  PDM_quaternion_t* qt
)
{
  double len = PDM_quaternion_norm(qt);
  qt->w  = qt->w  / len;
  qt->v[0] = qt->v[0] / len;
  qt->v[1] = qt->v[1] / len;
  qt->v[2] = qt->v[2] / len;
}

void
PDM_quaternion_compose
(
  const PDM_quaternion_t* qt_1,
  const PDM_quaternion_t* qt_2,
  PDM_quaternion_t* qt_out
)
{
  // buffering data in case where qt_1==qt_out or qt_2==qt_out
  double w,v1,v2,v3;
  w  = qt_1->w  * qt_2->w  - qt_1->v[0] * qt_2->v[0] - qt_1->v[1] * qt_2->v[1] - qt_1->v[2] * qt_2->v[2];
  v1 = qt_1->v[0] * qt_2->w  + qt_1->w  * qt_2->v[0] + qt_1->v[1] * qt_2->v[2] - qt_1->v[2] * qt_2->v[1];
  v2 = qt_1->w  * qt_2->v[1] - qt_1->v[0] * qt_2->v[2] + qt_1->v[1] * qt_2->w  + qt_1->v[2] * qt_2->v[0];
  v3 = qt_1->w  * qt_2->v[2] + qt_1->v[0] * qt_2->v[1] - qt_1->v[1] * qt_2->v[0] + qt_1->v[2] * qt_2->w ;
  qt_out->w  = w ;
  qt_out->v[0] = v1;
  qt_out->v[1] = v2;
  qt_out->v[2] = v3;
}

PDM_bool_t
PDM_quaternion_equal
(
  const PDM_quaternion_t* qt,
  const double w,
  const double v0,
  const double v1,
  const double v2,
  const double epsilon
)
{
  return (PDM_ABS(qt->w    - w ) < epsilon) & \
         (PDM_ABS(qt->v[0] - v0) < epsilon) &\
         (PDM_ABS(qt->v[1] - v1) < epsilon) &\
         (PDM_ABS(qt->v[2] - v2) < epsilon);
}

PDM_bool_t
PDM_quaternion_equal_quaternion
(
  const PDM_quaternion_t* qt_1,
  const PDM_quaternion_t* qt_2,
  const double epsilon
)
{
  return PDM_quaternion_equal(qt_1,qt_2->w,qt_2->v[0],qt_2->v[1],qt_2->v[2],epsilon);
}


/*----------------------------------------------------------------------------
 *  APPLICATION TO VECTORS
 *----------------------------------------------------------------------------*/

void
PDM_quaternion_rotate
(
  const PDM_quaternion_t* qt,
  const double* vector,
  const int n_samp,
  double* vector_out
)
{
  double ww = qt->w  * qt->w;
  double xx = qt->v[0] * qt->v[0];
  double yy = qt->v[1] * qt->v[1];
  double zz = qt->v[2] * qt->v[2];
  double wx = qt->w  * qt->v[0];
  double wy = qt->w  * qt->v[1];
  double wz = qt->w  * qt->v[2];
  double xy = qt->v[0] * qt->v[1];
  double xz = qt->v[0] * qt->v[2];
  double yz = qt->v[1] * qt->v[2];

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
  const PDM_quaternion_t* qt,
  const PDM_quaternion_t* qt_der,
  const double* vector,
  const int n_samp,
  double* vector_output_der
)
{
  double ww_der = 2*qt->w  * qt_der->w;
  double xx_der = 2*qt->v[0] * qt_der->v[0];
  double yy_der = 2*qt->v[1] * qt_der->v[1];
  double zz_der = 2*qt->v[2] * qt_der->v[2];
  double wx_der = qt_der->w  * qt->v[0] + qt->w  * qt_der->v[0];
  double wy_der = qt_der->w  * qt->v[1] + qt->w  * qt_der->v[1];
  double wz_der = qt_der->w  * qt->v[2] + qt->w  * qt_der->v[2];
  double xy_der = qt_der->v[0] * qt->v[1] + qt->v[0] * qt_der->v[1];
  double xz_der = qt_der->v[0] * qt->v[2] + qt->v[0] * qt_der->v[2];
  double yz_der = qt_der->v[1] * qt->v[2] + qt->v[1] * qt_der->v[2];

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
  PDM_quaternion_t qt;
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
  const double vector_1[3],
  const double vector_2[3],
  const double* t,
  const int n_samp,
  double* out
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_t qt_der;
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

void 
PDM_quaternion_from_two_vectors
(
  const double* vector_1,
  const double* vector_2,
  PDM_quaternion_t* qt_out
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
    qt_out->w = 1+dot;
    qt_out->v[0] = cross[0];
    qt_out->v[1] = cross[1];
    qt_out->v[2] = cross[2];
    PDM_quaternion_normalize(qt_out);
  }
}

void 
PDM_quaternion_from_axis_angle
(
  const double axis[3],
  const double angle,
  PDM_quaternion_t* qt_out
)
{ 
  qt_out->w = cos(0.5*angle);
  double c  = sin(0.5*angle);
  double axis_invnorm = 1./sqrt(PDM_DOT_PRODUCT(axis,axis));
  qt_out->v[0] = c * axis[0]*axis_invnorm;
  qt_out->v[1] = c * axis[1]*axis_invnorm;
  qt_out->v[2] = c * axis[2]*axis_invnorm;
}

void
PDM_quaternion_from_axis_angle_derivative(
  const double axis[3], 
  const double angle,
  const double axis_der[3],
  const double angle_der,
  PDM_quaternion_t* qt_out
)
{
  double c     =      sin(0.5*angle);
  qt_out->w    = -0.5*sin(0.5*angle)*angle_der;
  double c_der =  0.5*cos(0.5*angle)*angle_der;
  double ax_sq_norm = PDM_DOT_PRODUCT(axis,axis);
  double axis_invnorm     = 1./sqrt(ax_sq_norm);
  double axis_invnorm_der = -(axis[0]*axis_der[0]+axis[1]*axis_der[1]+axis[2]*axis_der[2])/(ax_sq_norm*sqrt(ax_sq_norm));
  qt_out->v[0] = c_der*axis[0]*axis_invnorm + c*axis_der[0]*axis_invnorm+c*axis[0]*axis_invnorm_der;
  qt_out->v[1] = c_der*axis[1]*axis_invnorm + c*axis_der[1]*axis_invnorm+c*axis[1]*axis_invnorm_der;
  qt_out->v[2] = c_der*axis[2]*axis_invnorm + c*axis_der[2]*axis_invnorm+c*axis[2]*axis_invnorm_der;
}

void
PDM_quaternion_from_axis_aligned_rotation
(
  const double angle,
  const int axis_ind,
  PDM_quaternion_t* qt_out
)
{
  qt_out->w  = cos(angle*0.5);
  qt_out->v[0] = 0.;
  qt_out->v[1] = 0.;
  qt_out->v[2] = 0.;
  qt_out->v[axis_ind] = sin(angle*0.5);
}

void
PDM_quaternion_from_x_rotation
(
  const double angle,
  PDM_quaternion_t* qt
)
{
  PDM_quaternion_from_axis_aligned_rotation(angle,0,qt);
}

void
PDM_quaternion_from_y_rotation
(
  const double angle,
  PDM_quaternion_t* qt
)
{
  PDM_quaternion_from_axis_aligned_rotation(angle,1,qt);
}

void
PDM_quaternion_from_z_rotation
(
  const double angle,
  PDM_quaternion_t* qt
)
{
  PDM_quaternion_from_axis_aligned_rotation(angle,2,qt);
}

void
PDM_quaternion_from_axis_aligned_symmetry
(
  const int axis_ind,
  PDM_quaternion_t* qt
)
{
  qt->w  = 0.;
  qt->v[0] = 0.;
  qt->v[1] = 0.;
  qt->v[2] = 0.;
  qt->v[axis_ind] = 1.;
}

void
PDM_quaternion_from_x_symmetry
(
  PDM_quaternion_t* qt
)
{
  PDM_quaternion_from_axis_aligned_symmetry(0,qt);
}

void
PDM_quaternion_from_y_symmetry
(
  PDM_quaternion_t* qt
)
{
  PDM_quaternion_from_axis_aligned_symmetry(1,qt);
}

void
PDM_quaternion_from_z_symmetry
(
  PDM_quaternion_t* qt
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
  PDM_quaternion_t* qt
)
{
  double angles[3] = {ang_x,ang_y,ang_z};
  PDM_quaternion_set_identity(qt);
  PDM_quaternion_t q_curr;
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
  const double *rotation_matrix,
  PDM_quaternion_t* qt_out
)
{
  double tr = rotation_matrix[3*0+0]+rotation_matrix[3*1+1]+rotation_matrix[3*2+2];
  double tmp;
  if (PDM_ABS(tr) > QUATERNION_EPS){
    tmp = sqrt(1+tr)*2.;
    qt_out->w    = .25*tmp;
    qt_out->v[0] = (rotation_matrix[3*2+1]-rotation_matrix[3*1+2])/tmp;
    qt_out->v[1] = (rotation_matrix[3*0+2]-rotation_matrix[3*2+0])/tmp;
    qt_out->v[2] = (rotation_matrix[3*1+0]-rotation_matrix[3*0+1])/tmp;
  }
  else if ((rotation_matrix[3*0+0]>rotation_matrix[3*1+1]) & (rotation_matrix[3*0+0]>rotation_matrix[3*2+2])){
    tmp = sqrt(1+rotation_matrix[3*0+0]-rotation_matrix[3*1+1]-rotation_matrix[3*2+2])*2.;
    qt_out->w    = (rotation_matrix[3*2+1]-rotation_matrix[3*1+2])/tmp;
    qt_out->v[0] = .25*tmp;
    qt_out->v[1] = (rotation_matrix[3*1+0]+rotation_matrix[3*0+1])/tmp;
    qt_out->v[2] = (rotation_matrix[3*2+0]+rotation_matrix[3*0+2])/tmp;
  }
  else if (rotation_matrix[3*1+1]>rotation_matrix[3*2+2]){
    tmp  = sqrt(1+rotation_matrix[3*1+1] - rotation_matrix[3*0+0] - rotation_matrix[3*2+2]) * 2;
    qt_out->w    = (rotation_matrix[3*0+2] - rotation_matrix[3*2+0]) / tmp;
    qt_out->v[0] = (rotation_matrix[3*0+1] + rotation_matrix[3*1+0]) / tmp; 
    qt_out->v[1] = 0.25 * tmp;
    qt_out->v[2] = (rotation_matrix[3*1+2] + rotation_matrix[3*2+1]) / tmp; 
  } else { 
    tmp = sqrt(1.0 + rotation_matrix[3*2+2] - rotation_matrix[3*0+0] - rotation_matrix[3*1+1]) * 2;
    qt_out->w    = (rotation_matrix[3*1+0] - rotation_matrix[3*0+1]) / tmp;
    qt_out->v[0] = (rotation_matrix[3*0+2] + rotation_matrix[3*2+0]) / tmp;
    qt_out->v[1] = (rotation_matrix[3*1+2] + rotation_matrix[3*2+1]) / tmp;
    qt_out->v[2] = 0.25 * tmp;
  }
}

void
PDM_quaternion_from_homogeneous_matrix
(
  const double *homogeneous_matrix,
  PDM_quaternion_t* qt_out
)
{
  double tr = homogeneous_matrix[4*0+0]+homogeneous_matrix[4*1+1]+homogeneous_matrix[4*2+2];
  double tmp;
  if (PDM_ABS(tr) > QUATERNION_EPS){
    tmp = sqrt(1+tr)*2.;
    qt_out->w    = .25*tmp;
    qt_out->v[0] = (homogeneous_matrix[4*2+1]-homogeneous_matrix[4*1+2])/tmp;
    qt_out->v[1] = (homogeneous_matrix[4*0+2]-homogeneous_matrix[4*2+0])/tmp;
    qt_out->v[2] = (homogeneous_matrix[4*1+0]-homogeneous_matrix[4*0+1])/tmp;
  }
  else if ((homogeneous_matrix[4*0+0]>homogeneous_matrix[4*1+1]) & (homogeneous_matrix[4*0+0]>homogeneous_matrix[4*2+2])){
    tmp = sqrt(1+homogeneous_matrix[4*0+0]-homogeneous_matrix[4*1+1]-homogeneous_matrix[4*2+2])*2.;
    qt_out->w    = (homogeneous_matrix[4*2+1]-homogeneous_matrix[4*1+2])/tmp;
    qt_out->v[0] = .25*tmp;
    qt_out->v[1] = (homogeneous_matrix[4*1+0]+homogeneous_matrix[4*0+1])/tmp;
    qt_out->v[2] = (homogeneous_matrix[4*2+0]+homogeneous_matrix[4*0+2])/tmp;
  }
  else if (homogeneous_matrix[4*1+1]>homogeneous_matrix[4*2+2]){
    tmp  = sqrt(1+homogeneous_matrix[4*1+1] - homogeneous_matrix[4*0+0] - homogeneous_matrix[4*2+2]) * 2;
    qt_out->w    = (homogeneous_matrix[4*0+2] - homogeneous_matrix[4*2+0]) / tmp;
    qt_out->v[0] = (homogeneous_matrix[4*0+1] + homogeneous_matrix[4*1+0]) / tmp; 
    qt_out->v[1] = 0.25 * tmp;
    qt_out->v[2] = (homogeneous_matrix[4*1+2] + homogeneous_matrix[4*2+1]) / tmp; 
  } else { 
    tmp = sqrt(1.0 + homogeneous_matrix[4*2+2] - homogeneous_matrix[4*0+0] - homogeneous_matrix[4*1+1]) * 2;
    qt_out->w    = (homogeneous_matrix[4*1+0] - homogeneous_matrix[4*0+1]) / tmp;
    qt_out->v[0] = (homogeneous_matrix[4*0+2] + homogeneous_matrix[4*2+0]) / tmp;
    qt_out->v[1] = (homogeneous_matrix[4*1+2] + homogeneous_matrix[4*2+1]) / tmp;
    qt_out->v[2] = 0.25 * tmp;
  }
}

/*----------------------------------------------------------------------------
 *  PDM_quaternion_t TO OTHER FORMATS
 *----------------------------------------------------------------------------*/

void
PDM_quaternion_to_axis_angle
(
  const PDM_quaternion_t* qt,
  double  axis[3],
  double* angle
)
{
  (*angle) = 2.0 * acos(qt->w);
  double divider = sqrt(1.0 - qt->w * qt->w);

  if(PDM_ABS(divider) > QUATERNION_EPS) {
    // Calculate the axis
    axis[0] = qt->v[0] / divider;
    axis[1] = qt->v[1] / divider;
    axis[2] = qt->v[2] / divider;
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
  const PDM_quaternion_t* qt,
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
    double quat[4] = {qt->v[0],qt->v[1],qt->v[2],qt->w};
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

    (*ang_x) = angle[order[0]];
    (*ang_y) = angle[order[1]];
    (*ang_z) = angle[order[2]];

}

void
PDM_quaternion_to_rotation_matrix
(
  const PDM_quaternion_t* qt, 
  double* rotation_matrix
)
{
  rotation_matrix[3*0+0] = 1 - 2 * qt->v[1] * qt->v[1] - 2 * qt->v[2] * qt->v[2];
  rotation_matrix[3*0+1] =     2 * qt->v[0] * qt->v[1] - 2 * qt->v[2] * qt->w ;
  rotation_matrix[3*0+2] =     2 * qt->v[0] * qt->v[2] + 2 * qt->v[1] * qt->w ;
  rotation_matrix[3*1+0] =     2 * qt->v[0] * qt->v[1] + 2 * qt->v[2] * qt->w ;
  rotation_matrix[3*1+1] = 1 - 2 * qt->v[2] * qt->v[2] - 2 * qt->v[0] * qt->v[0];
  rotation_matrix[3*1+2] =     2 * qt->v[1] * qt->v[2] - 2 * qt->v[0] * qt->w ;
  rotation_matrix[3*2+0] =     2 * qt->v[0] * qt->v[2] - 2 * qt->v[1] * qt->w ;
  rotation_matrix[3*2+1] =     2 * qt->v[1] * qt->v[2] + 2 * qt->v[0] * qt->w ;
  rotation_matrix[3*2+2] = 1 - 2 * qt->v[0] * qt->v[0] - 2 * qt->v[1] * qt->v[1];
}

void
PDM_quaternion_to_homogeneous_matrix
(
  const PDM_quaternion_t* qt,
  double* homogeneous_matrix
)
{
  homogeneous_matrix[4*0+0] = 1 - 2 * qt->v[1] * qt->v[1] - 2 * qt->v[2] * qt->v[2];
  homogeneous_matrix[4*0+1] =     2 * qt->v[0] * qt->v[1] - 2 * qt->v[2] * qt->w ;
  homogeneous_matrix[4*0+2] =     2 * qt->v[0] * qt->v[2] + 2 * qt->v[1] * qt->w ;
  homogeneous_matrix[4*0+3] = 0.;
  homogeneous_matrix[4*1+0] =     2 * qt->v[0] * qt->v[1] + 2 * qt->v[2] * qt->w ;
  homogeneous_matrix[4*1+1] = 1 - 2 * qt->v[2] * qt->v[2] - 2 * qt->v[0] * qt->v[0];
  homogeneous_matrix[4*1+2] =     2 * qt->v[1] * qt->v[2] - 2 * qt->v[0] * qt->w ;
  homogeneous_matrix[4*1+3] = 0.;
  homogeneous_matrix[4*2+0] =     2 * qt->v[0] * qt->v[2] - 2 * qt->v[1] * qt->w ;
  homogeneous_matrix[4*2+1] =     2 * qt->v[1] * qt->v[2] + 2 * qt->v[0] * qt->w ;
  homogeneous_matrix[4*2+2] = 1 - 2 * qt->v[0] * qt->v[0] - 2 * qt->v[1] * qt->v[1];
  homogeneous_matrix[4*2+3] = 0.;
  homogeneous_matrix[4*3+0] = 0.;
  homogeneous_matrix[4*3+1] = 0.;
  homogeneous_matrix[4*3+2] = 0.;
  homogeneous_matrix[4*3+3] = 1.;
}

/*----------------------------------------------------------------------------
 *  PDM_quaternion_t INPUT/OUTPUT FUNCTIONS
 *----------------------------------------------------------------------------*/

void 
PDM_quaternion_axis_angle_to_euler_angles
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
  PDM_quaternion_t qt;
  PDM_quaternion_from_axis_angle(axis,angle,&qt);
  PDM_quaternion_to_euler_angles(&qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void 
PDM_quaternion_axis_angle_to_rotation_matrix
(
  const double axis[3],
  const double angle,
  double *rotation_matrix
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_from_axis_angle(axis,angle,&qt);
  PDM_quaternion_to_rotation_matrix(&qt,rotation_matrix);
}

void 
PDM_quaternion_axis_angle_to_homogeneous_matrix
(
  const double axis[3],
  const double angle,
  double *homogeneous_matrix
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_from_axis_angle(axis,angle,&qt);
  PDM_quaternion_to_homogeneous_matrix(&qt,homogeneous_matrix);
}

void
PDM_quaternion_euler_angles_to_axis_angle
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
  PDM_quaternion_t qt;
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&qt);
  PDM_quaternion_to_axis_angle(&qt,axis,angle);
}

void
PDM_quaternion_euler_angles_to_euler_angles
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
  PDM_quaternion_t qt;
  PDM_quaternion_from_euler_angles(input_ang_x,input_ang_y,input_ang_z,input_order,input_intrinsic,&qt);
  PDM_quaternion_to_euler_angles(&qt,output_order,output_intrinsic,output_ang_x,output_ang_y,output_ang_z);
}

void
PDM_quaternion_euler_angles_to_rotation_matrix
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  PDM_bool_t intrinsic,
  double* rotation_matrix
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&qt);
  PDM_quaternion_to_rotation_matrix(&qt,rotation_matrix);
}

void
PDM_quaternion_euler_angles_to_homogeneous_matrix
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  PDM_bool_t intrinsic,
  double* homogeneous_matrix
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_from_euler_angles(ang_x,ang_y,ang_z,order,intrinsic,&qt);
  PDM_quaternion_to_homogeneous_matrix(&qt,homogeneous_matrix);
}


void 
PDM_quaternion_rotation_matrix_to_axis_angle
(
  const double* rotation_matrix,
  double axis[3],
  double* angle
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_from_rotation_matrix(rotation_matrix,&qt);
  PDM_quaternion_to_axis_angle(&qt,axis,angle);

}

void 
PDM_quaternion_rotation_matrix_to_euler_angles
(
  const double* rotation_matrix,
  const int order[3],
  const PDM_bool_t intrinsic,
  double* ang_x,
  double* ang_y,
  double* ang_z
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_from_rotation_matrix(rotation_matrix,&qt);
  PDM_quaternion_to_euler_angles(&qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void 
PDM_quaternion_rotation_matrix_to_homogeneous_matrix
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
PDM_quaternion_homogeneous_matrix_to_axis_angle
(
  const double* homogeneous_matrix,
  double axis[3],
  double* angle
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_from_homogeneous_matrix(homogeneous_matrix,&qt);
  PDM_quaternion_to_axis_angle(&qt,axis,angle);
}

void 
PDM_quaternion_homogeneous_matrix_to_euler_angles
(
  const double* homogeneous_matrix,
  const int order[3],
  const PDM_bool_t intrinsic,
  double* ang_x,
  double* ang_y,
  double* ang_z
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_from_homogeneous_matrix(homogeneous_matrix,&qt);
  PDM_quaternion_to_euler_angles(&qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void 
PDM_quaternion_homogeneous_matrix_to_rotation_matrix
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
PDM_quaternion_two_vectors_to_axis_angle
(
  const double vector_1[3],
  const double vector_2[3],
  double axis[3],
  double* angle
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,&qt);
  PDM_quaternion_to_axis_angle(&qt,axis,angle);
}

void 
PDM_quaternion_two_vectors_to_euler_angles
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
  PDM_quaternion_t qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,&qt);
  PDM_quaternion_to_euler_angles(&qt,order,intrinsic,ang_x,ang_y,ang_z);
}

void 
PDM_quaternion_two_vectors_to_rotation_matrix
(
  const double vector_1[3],
  const double vector_2[3],
  double *rotation_matrix
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,&qt);
  PDM_quaternion_to_rotation_matrix(&qt,rotation_matrix);
}

void 
PDM_quaternion_two_vectors_to_homogeneous_matrix
(
  const double vector_1[3],
  const double vector_2[3],
  double *homogeneous_matrix
)
{
  PDM_quaternion_t qt;
  PDM_quaternion_from_two_vectors(vector_1,vector_2,&qt);
  PDM_quaternion_to_homogeneous_matrix(&qt,homogeneous_matrix);
}

/*----------------------------------------------------------------------------
 *  PDM_quaternion_t COMPOSITE FUNCTIONS
 *----------------------------------------------------------------------------*/

#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
char BlasNoTrans = 'N';
char BlasTrans = 'T';
extern void dgemm_(char  	*TransA,
                  char  	*TransB,
                  int   	*M,
                  int   	*N,
                  int   	*K,
                  double * 	alpha,
                  double *  	A,
                  int  	*lda,
                  double *  	B,
                  int  	*ldb,
                  double * 	beta,
                  double *  	C,
                  int  	*ldc);
#endif

void 
PDM_quaternion_identity_to_homogeneous_matrix
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


void
PDM_quaternion_multiply_homogeneous_matrices
(
  const double A[16],
  const double B[16],
  double       C[16]
)
{
#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
  int ldABC = 4;
  int MNK = 4;
  double alpha = 1.;
  double beta = 0.;
  dgemm_(&BlasNoTrans,&BlasNoTrans,&MNK,&MNK,&MNK,&alpha,A,&ldABC,
    B,&ldABC,&beta,C,&ldABC);
#else
  printf("Error : LAPACK or MKL are mandatory, recompile with them. \n");
  abort();
#endif
}

void 
PDM_quaternion_apply_rotation_matrix
(
  const double rotation_matrix[9],
  const double* vector,
  const int n_samp,
  double* vector_out
)
{
#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
  int ldB = 3;
  int ldA = n_samp;
  double alpha = 1.;
  double beta = 0.;
  int NK = 3;
  // dgemm_(&BlasNoTrans,&BlasTrans,&n_samp,&NK,&NK,&alpha,vector,&ldA,
  //       rotation_matrix,&ldB,&beta,vector_out,&ldA);
  dgemm_(&BlasNoTrans,&BlasNoTrans,&NK,&n_samp,&NK,&alpha,rotation_matrix,&ldB,
        vector,&ldB,&beta,vector_out,&ldB);
#else
  printf("Error : LAPACK or MKL are mandatory, recompile with them. \n");
  abort();
#endif
}

void
PDM_quaternion_apply_translation
(
  const double translation_vector[3],
  const double* vector,
  const int n_samp,
  double* vector_out
)
{
  for (int i = 0; i < n_samp; i++) {
    vector_out[3*i+0] = vector[3*i+0] + translation_vector[0];
    vector_out[3*i+1] = vector[3*i+1] + translation_vector[1];
    vector_out[3*i+2] = vector[3*i+2] + translation_vector[2];
  }
}

void 
PDM_quaternion_apply_homogeneous_matrix
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
  PDM_quaternion_apply_rotation_matrix(rotation_matrix,vector,n_samp,vector_out);

  // applying the translation
  for (int i = 0; i < n_samp; i++) {
    vector_out[3*i+0] = homogeneous_matrix[3];
    vector_out[3*i+1] = homogeneous_matrix[7];
    vector_out[3*i+2] = homogeneous_matrix[11];
  }
}

void
PDM_quaternion_compose_homogeneous_matrices
(
  const double** homogeneous_matrices,
  const int n_matrices,
  double output_matrix[16]
)
{
  PDM_quaternion_identity_to_homogeneous_matrix(output_matrix);
  for (int i = 0; i < n_matrices; i++) {
    PDM_quaternion_multiply_homogeneous_matrices(homogeneous_matrices[i],output_matrix,output_matrix);
  }
}

void
PDM_quaternion_translation_to_homogeneous_matrix
(
  const double translation_vector[3],
  PDM_bool_t reverse,
  double homogeneous_matrix[16]
)
{
  PDM_quaternion_identity_to_homogeneous_matrix(homogeneous_matrix);
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

void 
PDM_quaternion_apply_euler_angles_and_rotation_center
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  PDM_bool_t intrinsic,
  const double rotation_center[3],
  PDM_bool_t reverse,
  const double* vector,
  const int n_samp,
  double* vector_out
)
{
  // building the homogeneous matrix corresponding to 
  double homogeneous_matrix[16];
  double tmp_matrix[16];
  // translation of -rotation_center
  PDM_quaternion_translation_to_homogeneous_matrix(rotation_center,PDM_TRUE,homogeneous_matrix);
  // rotation
  PDM_quaternion_euler_angles_to_homogeneous_matrix(ang_x,ang_y,ang_z,order,intrinsic,tmp_matrix);
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
  PDM_quaternion_multiply_homogeneous_matrices(homogeneous_matrix,tmp_matrix,homogeneous_matrix);
  
  // // translation of rotation_center
  // PDM_quaternion_translation_to_homogeneous_matrix(rotation_center,PDM_FALSE,tmp_matrix);
  // PDM_quaternion_multiply_homogeneous_matrices(homogeneous_matrix,tmp_matrix,homogeneous_matrix);

  // // applying the homogeneous matrix to the coordinate vector
  // PDM_quaternion_apply_homogeneous_matrix(homogeneous_matrix,vector,n_samp,vector_out);

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
