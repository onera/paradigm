#ifndef __PDM_QUATERNION_H__
#define __PDM_QUATERNION_H__

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_quaternion_t PDM_quaternion_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 *
 * \brief Creates a quaternion
 *
 * \param [in]   w     Real (scalar) part
 * \param [in]   v1    First component of the imaginary (vectorial) part
 * \param [in]   v2    Second component of the imaginary (vectorial) part
 * \param [in]   v3    Third component of the imaginary (vectorial) part
 * \return The corresponding quaternion
 *
 */

PDM_quaternion_t*
PDM_quaternion_create
(
  const double w,
  const double v1,
  const double v2,
  const double v3
);

void
PDM_quaternion_free
(
  PDM_quaternion_t* qt
);

/**
 *
 * \brief Sets a quaternion to identity (1.,(0.,0.,0.))
 *
 * \param [in]   qt    Quaternion
 *
 */


void
PDM_quaternion_set_identity
(
  PDM_quaternion_t* qt
);

/**
 *
 * \brief Prints a quaternion
 *
 * \param [in]   qt    Quaternion
 *
 */

void
PDM_quaternion_print
(
  const PDM_quaternion_t* qt
);

/**
 *
 * \brief Computes the conjugate of a quaternion
 *
 * \param [in]   qt        Quaternion
 * \param [in]   qt_out    Conjugate quaternion
 *
 */

void
PDM_quaternion_conjugate
(
  const PDM_quaternion_t* qt,
  PDM_quaternion_t* qt_out
);


/**
 *
 * \brief Computes the norm of a quaternion (sqrt(w*w+v0*v0+v1*v1+v2*v2))
 *
 * \param [in]   qt        Quaternion
 * \return The L2-norm of the quaternion (double)
 *
 */

double
PDM_quaternion_norm
(
  const PDM_quaternion_t* qt
);

/**
 *
 * \brief Noramlizes a quaternion (q/sqrt(w*w+v0*v0+v1*v1+v2*v2))
 *
 * \param [in]   qt        Quaternion
 *
 */

void
PDM_quaternion_normalize
(
  PDM_quaternion_t* qt
);

/**
 *
 * \brief Composes 2 quaternions (quaternion multiplication)
 * Formula from http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/arithmetic/index.htm
 *        a*e - b*f - c*g - d*h
 *   + i (b*e + a*f + c*h - d*g)
 *   + j (a*g - b*h + c*e + d*f)
 *   + k (a*h + b*g - c*f + d*e)
 * 
 * /!\ Quaternion multiplication is not commutative !
 *
 * \param [in]   qt_1        First quaternion
 * \param [in]   qt_2        First quaternion
 * \return The result of qt_1*qt_2 
 *
 */

void
PDM_quaternion_compose
(
  const PDM_quaternion_t* qt_1,
  const PDM_quaternion_t* qt_2,
  PDM_quaternion_t* qt_out
);

/**
 *
 * \brief Tests the equality of a quaternion with provided components
 *
 * \param [in]   qt       Quaternion
 * \param [in]   w        Scalar part
 * \param [in]   v0       First component of the vector part
 * \param [in]   v1       Second component of the vector part
 * \param [in]   v2       Third component of the vector part
 * \return PDM_True in case of equality, PDM_False otherwise 
 *
 */

PDM_bool_t
PDM_quaternion_equal
(
  const PDM_quaternion_t* qt,
  const double w,
  const double v0,
  const double v1,
  const double v2,
  const double epsilon
);

/**
 *
 * \brief Tests the equality of 2 quaternions
 *
 * \param [in]   qt_1     First quaternion
 * \param [in]   qt_2     Second quaternion
 * \return PDM_True in case of equality, PDM_False otherwise 
 *
 */

PDM_bool_t
PDM_quaternion_equal_quaternion
(
  const PDM_quaternion_t* qt_1,
  const PDM_quaternion_t* qt_2,
  const double epsilon
);

/**
 *
 * \brief Applies a quaternion to (a) flatten 3D vector(s)
 *
 * \param [in]   qt           Quaternion
 * \param [in]   vector       The flatten vector array of length 3*n_samp
 * \param [in]   n_samp       The number of vector samples
 * \param [in]   vector_out   The rotated flatten vector array of length 3*n_samp
 *
 */

void
PDM_quaternion_rotate
(
  const PDM_quaternion_t* qt,
  const double* vector,
  const int n_samp,
  double* vector_out
);

/**
 *
 * \brief Computes the derivative of the rotation by qt in the direction qt_der
 *
 * \param [in]   qt                  Quaternion
 * \param [in]   qt_der              Differential of the quaternion
 * \param [in]   vector              The flatten vector array of length 3*n_samp
 * \param [in]   n_samp              The number of vector samples
 * \param [in]   vector_output_der   The derivative of the rotated flatten vector array of length 3*n_samp wrt to qt_der 
 *
 */

void
PDM_quaternion_rotate_derivative
(
  const PDM_quaternion_t* qt,
  const PDM_quaternion_t* qt_der,
  const double* vector,
  const int n_samp,
  double* vector_output_der
);

/**
 *
 * \brief Computes the spherical linear interpolation between two unit 3D vectors 
 *
 * \param [in]   vector_1            First vector
 * \param [in]   vector_2            Second vector
 * \param [in]   t                   The flatten array of t values of length n_samp (extrapolation is possible)
 * \param [in]   n_samp              The number of samples
 * \param [in]   out                 The result of the slerp of length 3*n_samp
 *
 */

void
PDM_quaternion_slerp_from_two_vectors
(
  const double vector_1[3],
  const double vector_2[3],
  const double* t,
  const int n_samp,
  double* out
);

/**
 *
 * \brief Computes the derivative of the spherical linear interpolation between two unit vectors with respect to t
 *
 * \param [in]   vector_1            First vector
 * \param [in]   vector_2            Second vector
 * \param [in]   t                   The flatten array of t values of length n_samp (extrapolation is possible)
 * \param [in]   n_samp              The number of samples
 * \param [in]   out                 The derivative of the slerp of length 3*n_samp
 *
 */

void
PDM_quaternion_slerp_from_two_vectors_derivative
(
  const double vector_1[3],
  const double vector_2[3],
  const double* t,
  const int n_samples,
  double* out
);

/**
 *
 * \brief Computes the quaternion corresponding to the rotation from the first 3D unit vector to the second one
 * Formula from https://stackoverflow.com/questions/1171849/finding-PDM_quaternion_t-representing-the-rotation-from-one-vector-to-another
 *
 * \param [in]   vector_1            First vector
 * \param [in]   vector_2            Second vector
 * \param [in]   qt_out              The output quaternion
 *
 */

void 
PDM_quaternion_from_two_vectors
(
  const double vector_1[3],
  const double vector_2[3],
  PDM_quaternion_t* qt_out
);

/**
 *
 * \brief Computes the quaternion corresponding to a rotation around the provided axis
 * Formula from http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/
 * \param [in]   axis            Rotation axis (3D (unit) vector)
 * \param [in]   angle           Rotation angle (in radians)
 * \param [in]   qt_out          The output quaternion
 *
 */

void 
PDM_quaternion_from_axis_angle
(
  const double axis[3],
  const double angle,
  PDM_quaternion_t* qt_out
);

/**
 *
 * \brief Computes the derivative of the quaternion corresponding to a rotation around the provided axis with respect to the angle and the axis
 * Formula from http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/
 *
 * \param [in]   axis            Rotation axis (3D (unit) vector)
 * \param [in]   angle           Rotation angle (in radians)
 * \param [in]   axis_der        Differential of the rotation axis (3D (unit) vector)
 * \param [in]   angle_der       Differential of the rotation angle (in radians)
 * \param [in]   qt_out          The output quaternion
 *
 */

void
PDM_quaternion_from_axis_angle_derivative(
  const double axis[3], 
  const double angle,
  const double axis_der[3],
  const double angle_der,
  PDM_quaternion_t* qt_out
);

/**
 *
 * \brief Computes the quaternion corresponding to a rotation around the x, y or z axis
 *
 * \param [in]   angle           Rotation angle (in radians)
 * \param [in]   axis_ind        Index of the axis (0 -> x, 1 -> y, 2 -> z)
 * \param [in]   qt_out          The output quaternion
 *
 */

void
PDM_quaternion_from_axis_aligned_rotation
(
  const double angle,
  const int axis_ind,
  PDM_quaternion_t* qt_out
);

/**
 *
 * \brief Computes the quaternion corresponding to a rotation around the x axis
 *
 * \param [in]   angle           Rotation angle (in radians)
 * \param [in]   qt_out          The output quaternion
 *
 */

void
PDM_quaternion_from_x_rotation
(
  const double angle,
  PDM_quaternion_t* qt
);

/**
 *
 * \brief Computes the quaternion corresponding to a rotation around the y axis
 *
 * \param [in]   angle           Rotation angle (in radians)
 * \param [in]   qt_out          The output quaternion
 *
 */

void
PDM_quaternion_from_y_rotation
(
  const double angle,
  PDM_quaternion_t* qt
);

/**
 *
 * \brief Computes the quaternion corresponding to a rotation around the z axis
 *
 * \param [in]   angle           Rotation angle (in radians)
 * \param [in]   qt_out          The output quaternion
 *
 */

void
PDM_quaternion_from_z_rotation
(
  const double angle,
  PDM_quaternion_t* qt
);

/**
 *
 * \brief Computes the quaternion corresponding to an axial symmetry around the x, y or z axis
 *
 * \param [in]   axis_ind        Index of the axis (0 -> x, 1 -> y, 2 -> z)
 * \param [in]   qt_out          The output quaternion
 *
 */

void
PDM_quaternion_from_axis_aligned_symmetry
(
  const int axis_ind,
  PDM_quaternion_t* qt
);

/**
 *
 * \brief Computes the quaternion corresponding to an axial symmetry around the x axis
 *
 * \param [in]   qt_out          The output quaternion
 *
 */

void
PDM_quaternion_from_x_symmetry
(
  PDM_quaternion_t* qt
);

/**
 *
 * \brief Computes the quaternion corresponding to an axial symmetry around the y axis
 *
 * \param [in]   qt_out          The output quaternion
 *
 */

void
PDM_quaternion_from_y_symmetry
(
  PDM_quaternion_t* qt
);

/**
 *
 * \brief Computes the quaternion corresponding to an axial symmetry around the z axis
 *
 * \param [in]   qt_out          The output quaternion
 *
 */

void
PDM_quaternion_from_z_symmetry
(
  PDM_quaternion_t* qt
);

/**
 *
 * \brief Computes the quaternion corresponding to a rotation corresponding to euler angles
 *
 * Standard practice for aero is:
 * - ang_x = gamma (roll), ang_y = alpha (pitch, aoa), ang_z = beta (yaw)
 * - order = {2,1,0} (Z rotation ,then Y then X)
 * - intrinsic = PDM_True
 *
 * \param [in]   ang_x       Rotation angle around the x-axis
 * \param [in]   ang_y       Rotation angle around the y-axis
 * \param [in]   ang_z       Rotation angle around the z-axis
 * \param [in]   order       Order of rotations to apply 
 * \param [in]   intrinsic   Axis conventions (https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations)
 * \param [in]   qt          The output quaternion
 *
 */

void 
PDM_quaternion_from_euler_angles
(
  const double ang_x,
  const double ang_y,
  const double ang_z,
  const int order[3],
  PDM_bool_t intrinsic,
  PDM_quaternion_t* qt
);

/**
 *
 * \brief Computes the quaternion corresponding to a rotation matrix
 * Formula from http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
 *
 * \param [in]   rotation_matrix 3-by-3 rotation matrix
 * \param [in]   qt_out          The output quaternion
 *
 */

void
PDM_quaternion_from_rotation_matrix
(
  const double **rotation_matrix,
  PDM_quaternion_t* qt_out
);

/**
 *
 * \brief Computes the quaternion corresponding to a rotation matrix in homogeneous coordinates (4-by-4)
 *
 * \param [in]   homogeneous_matrix 4-by-4 homogeneous rotation matrix
 * \param [in]   qt_out          The output quaternion
 *
 */

void
PDM_quaternion_from_homogeneous_matrix
(
  const double **homogeneous_matrix,
  PDM_quaternion_t* qt_out
);

/**
 *
 * \brief Computes the axis and angle corresponding to the provided quaternion
 * Formula from http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/
 *
 * \param [in]   qt     Quaternion
 * \param [in]   axis   Rotation 3D-axis
 * \param [in]   angle  Rotation angle (in radians)
 *
 */


void
PDM_quaternion_to_axis_angle
(
  const PDM_quaternion_t* qt,
  double  axis[3],
  double* angle
);

/**
 *
 * \brief Computes the Euler angles corresponding to the provided quaternion
 * (shamelessly) stolen from https://github.com/scipy/scipy/blob/main/scipy/spatial/transform/_rotation.pyx#L359
 *
 * \param [in]   qt         Quaternion
 * \param [in]   order      Order of rotations
 * \param [in]   intrinsic  Axis conventions (https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations)
 * \param [in]   ang_x      Rotation angle around the x axis (in radians)
 * \param [in]   ang_y      Rotation angle around the y axis (in radians)
 * \param [in]   ang_z      Rotation angle around the z axis (in radians)
 *
 */

void
PDM_quaternion_to_euler_angles
(
  const PDM_quaternion_t* qt,
  const int order[3],
  const PDM_bool_t intrinsic,
  double* ang_x,
  double* ang_y,
  double* ang_z
);

/**
 *
 * \brief Computes the 3-by-3 rotation matrix corresponding to the provided quaternion
 * Formula from http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm
 *
 * \param [in]   qt                Quaternion
 * \param [in]   rotation_matrix   3-by-3 rotation matrix
 *
 */

void
PDM_quaternion_to_rotation_matrix
(
  const PDM_quaternion_t* qt, 
  double** rotation_matrix
);

/**
 *
 * \brief Computes the 3-by-3 rotation matrix corresponding to the provided quaternion
 *
 * \param [in]   qt                   Quaternion
 * \param [in]   homogeneous_matrix   4-by-4 homogeneous rotation matrix
 *
 */
void
PDM_quaternion_to_homogeneous_matrix
(
  const PDM_quaternion_t* qt,
  double** homogeneous_matrix
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_QUATERNION_H__ */
