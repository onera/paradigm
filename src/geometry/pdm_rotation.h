#ifndef __PDM_ROTATION_H__
#define __PDM_ROTATION_H__

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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/



/*----------------------------------------------------------------------------
 *  PDM_rotation_t INPUT/OUTPUT FUNCTIONS
 *----------------------------------------------------------------------------*/


// Axis angle to other formats ---

/**
 *
 * \brief Converts a rotation expressed as axis-angle to euler angles
 *
 * \param [in]   axis        Rotation axis (3D (unit) vector)
 * \param [in]   angle       Rotation angle (in radians)
 * \param [in]   reverse     If True, encodes the reverse transformation
 * \param [in]   order       Order of rotations to apply
 * \param [in]   intrinsic   <a href="https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations">Axis conventions</a>
 * \param [out]  ang_x       Rotation angle around the x-axis
 * \param [out]  ang_y       Rotation angle around the y-axis
 * \param [out]  ang_z       Rotation angle around the z-axis
 *
 */


void
PDM_rotation_axis_angle_to_euler_angles
(
  const double      axis[3],
  const double      angle,
  const PDM_bool_t  reverse,
  const int         order[3],
  const PDM_bool_t  intrinsic,
        double     *ang_x,
        double     *ang_y,
        double     *ang_z
);

/**
 *
 * \brief Converts a rotation expressed as axis-angle to a 3-by-3 rotation matrix
 *
 * \param [in]   axis            Rotation axis (3D (unit) vector)
 * \param [in]   angle           Rotation angle (in radians)
 * \param [in]   reverse         If True, encodes the reverse transformation
 * \param [out]  rotation_matrix 3-by-3 rotation matrix
 *
 */

void
PDM_rotation_axis_angle_to_rotation_matrix
(
  const double      axis[3],
  const double      angle,
  const PDM_bool_t  reverse,
        double     *rotation_matrix
);

/**
 *
 * \brief Converts a rotation expressed as axis-angle to a 4-by-4 homogeneous rotation matrix
 *
 * \param [in]   axis               Rotation axis (3D (unit) vector)
 * \param [in]   angle              Rotation angle (in radians)
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  homogeneous_matrix 4-by-4 rotation matrix
 *
 */

void
PDM_rotation_axis_angle_to_homogeneous_matrix
(
  const double      axis[3],
  const double      angle,
  const PDM_bool_t  reverse,
        double     *homogeneous_matrix
);

/**
 *
 * \brief Converts a rotation expressed as axis-angle and a rotation center
 * to a 4-by-4 homogeneous rotation matrix
 *
 * \param [in]   axis               Rotation axis (3D (unit) vector)
 * \param [in]   angle              Rotation angle (in radians)
 * \param [in]   rotation_center    3D rotation center
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  homogeneous_matrix 4-by-4 rotation matrix
 *
 */

void
PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix
(
  const double      axis[3],
  const double      angle,
  const double      rotation_center[3],
  const PDM_bool_t  reverse,
        double     *homogeneous_matrix
);

/**
 *
 * \brief Converts a rotation expressed as euler angles to an axis-angle rotation
 *
 * \param [in]   ang_x       Rotation angle around the x-axis
 * \param [in]   ang_y       Rotation angle around the y-axis
 * \param [in]   ang_z       Rotation angle around the z-axis
 * \param [in]   order       Order of rotations to apply
 * \param [in]   intrinsic   <a href="https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations">Axis conventions</a>
 * \param [in]   reverse     If True, encodes the reverse transformation
 * \param [out]  axis        Rotation axis (3D (unit) vector)
 * \param [out]  angle       Rotation angle (in radians)
 *
 */

void
PDM_rotation_euler_angles_to_axis_angle
(
  const double      ang_x,
  const double      ang_y,
  const double      ang_z,
  const int         order[3],
  const PDM_bool_t  intrinsic,
  const PDM_bool_t  reverse,
        double      axis[3],
        double     *angle
);

// Euler angles to other formats ---

/**
 *
 * \brief Converts a rotation expressed as euler angles to another euler angles representation
 *
 * \param [in]   input_ang_x       Input rotation angle around the x-axis
 * \param [in]   input_ang_y       Input rotation angle around the y-axis
 * \param [in]   input_ang_z       Input rotation angle around the z-axis
 * \param [in]   input_order       Input order of rotations to apply
 * \param [in]   input_intrinsic   Input <a href="https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations">axis conventions</a>
 * \param [in]   reverse           If True, encodes the reverse transformation
 * \param [in]   output_order      Output order of rotations to apply
 * \param [in]   output_intrinsic  Output <a href="https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations">axis conventions</a>
 * \param [out]  output_ang_x      Output rotation angle around the x-axis
 * \param [out]  output_ang_y      Output rotation angle around the y-axis
 * \param [out]  output_ang_z      Output rotation angle around the z-axis
 *
 */

void
PDM_rotation_euler_angles_to_euler_angles
(
  const double      input_ang_x,
  const double      input_ang_y,
  const double      input_ang_z,
  const int         input_order[3],
  const PDM_bool_t  input_intrinsic,
  const PDM_bool_t  reverse,
  const int         output_order[3],
  const PDM_bool_t  output_intrinsic,
        double     *output_ang_x,
        double     *output_ang_y,
        double     *output_ang_z
);

/**
 *
 * \brief Converts a rotation expressed as euler angles to a 3-by-3 rotation matrix
 *
 * \param [in]   ang_x            Rotation angle around the x-axis
 * \param [in]   ang_y            Rotation angle around the y-axis
 * \param [in]   ang_z            Rotation angle around the z-axis
 * \param [in]   order            Order of rotations to apply
 * \param [in]   intrinsic        <a href="https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations">Axis conventions</a>
 * \param [in]   reverse          If True, encodes the reverse transformation
 * \param [out]  rotation_matrix  3-by-3 rotation matrix
 *
 */

void
PDM_rotation_euler_angles_to_rotation_matrix
(
  const double      ang_x,
  const double      ang_y,
  const double      ang_z,
  const int         order[3],
  const PDM_bool_t  intrinsic,
  const PDM_bool_t  reverse,
        double     *rotation_matrix
);

/**
 *
 * \brief Converts a rotation expressed as euler angles to a 4-by-4 homogeneous rotation matrix
 *
 * \param [in]   ang_x              Rotation angle around the x-axis
 * \param [in]   ang_y              Rotation angle around the y-axis
 * \param [in]   ang_z              Rotation angle around the z-axis
 * \param [in]   order              Order of rotations to apply
 * \param [in]   intrinsic          <a href="https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations">Axis conventions</a>
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  homogeneous_matrix 4-by-4 rotation matrix
 *
 */

void
PDM_rotation_euler_angles_to_homogeneous_matrix
(
  const double      ang_x,
  const double      ang_y,
  const double      ang_z,
  const int         order[3],
  const PDM_bool_t  intrinsic,
  const PDM_bool_t  reverse,
        double     *homogeneous_matrix
);

/**
 *
 * \brief Converts a rotation expressed as euler angles and a rotation center
 * to a 4-by-4 homogeneous rotation matrix
 *
 * \param [in]   ang_x              Rotation angle around the x-axis
 * \param [in]   ang_y              Rotation angle around the y-axis
 * \param [in]   ang_z              Rotation angle around the z-axis
 * \param [in]   order              Order of rotations to apply
 * \param [in]   intrinsic          <a href="https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations">Axis conventions</a>
 * \param [in]   rotation_center    Rotation center
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  homogeneous_matrix 4-by-4 rotation matrix
 *
 */

void
PDM_rotation_euler_angles_and_rotation_center_to_homogeneous_matrix
(
  const double      ang_x,
  const double      ang_y,
  const double      ang_z,
  const int         order[3],
  const PDM_bool_t  intrinsic,
  const double      rotation_center[3],
  const PDM_bool_t  reverse,
        double     *homogeneous_matrix
);

/**
 *
 * \brief Converts the info of a CGNS Periodic_t node to a 4-by-4 homogeneous rotation matrix
 *
 * \param [in]   rotation_center
 * \param [in]   rotation_angle
 * \param [in]   translation
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  homogeneous_matrix 4-by-4 rotation matrix
 *
 */

void
PDM_rotation_periodic_t_info_to_homogeneous_matrix
(
  const double      rotation_center[3],
  const double      rotation_angle[3],
  const double      translation[3],
  const PDM_bool_t  reverse,
        double     *homogeneous_matrix
);

// Rotation matrix to other formats ---

/**
 *
 * \brief Converts a rotation expressed as a 3-by-3 rotation matrix to an axis-angle representation
 *
 * \param [in]   rotation_matrix 3-by-3 rotation matrix
 * \param [in]   reverse         If True, encodes the reverse transformation
 * \param [out]  axis            Rotation axis (3D (unit) vector)
 * \param [out]  angle           Rotation angle (in radians)
 *
 */

void
PDM_rotation_rotation_matrix_to_axis_angle
(
  const double     *rotation_matrix,
  const PDM_bool_t  reverse,
        double      axis[3],
        double     *angle
);

/**
 *
 * \brief Converts a rotation expressed as a 3-by-3 rotation matrix to euler angles
 *
 * \param [in]   rotation_matrix  3-by-3 rotation matrix
 * \param [in]   reverse          If True, encodes the reverse transformation
 * \param [in]   order            Order of rotations to apply
 * \param [in]   intrinsic        <a href="https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations">Axis conventions</a>
 * \param [out]  ang_x            Rotation angle around the x-axis
 * \param [out]  ang_y            Rotation angle around the y-axis
 * \param [out]  ang_z            Rotation angle around the z-axis
 *
 */

void
PDM_rotation_rotation_matrix_to_euler_angles
(
  const double     *rotation_matrix,
  const PDM_bool_t  reverse,
  const int         order[3],
  const PDM_bool_t  intrinsic,
        double     *ang_x,
        double     *ang_y,
        double     *ang_z
);

/**
 *
 * \brief Converts a rotation expressed as a 3-by-3 rotation matrix to a 4-by-4 homogeneous matrix
 *
 * \param [in]   rotation_matrix    3-by-3 rotation matrix
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  homogeneous_matrix 4-by-4 homogeneous rotation matrix
 *
 */

void
PDM_rotation_rotation_matrix_to_homogeneous_matrix
(
  const double     *rotation_matrix,
  const PDM_bool_t  reverse,
        double     *homogeneous_matrix
);

/**
 *
 * \brief Converts a rotation expressed as a 3-by-3 rotation matrix and a rotation center
 * to a 4-by-4 homogeneous matrix
 *
 * \param [in]   rotation_matrix    3-by-3 rotation matrix
 * \param [in]   rotation_center    Rotation center
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  homogeneous_matrix 4-by-4 homogeneous rotation matrix
 *
 */

void
PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix
(
  const double     *rotation_matrix,
  const double      rotation_center[3],
  const PDM_bool_t  reverse,
        double     *homogeneous_matrix
);

// Homogeneous matrix to other formats ---

/**
 *
 * \brief Converts a rotation expressed as a 4-by-4 homogeneous rotation matrix to an axis-angle representation
 *
 * \param [in]   homogeneous_matrix 4-by-4 homogeneous rotation matrix
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  axis               Rotation axis (3D (unit) vector)
 * \param [out]  angle              Rotation angle (in radians)
 *
 */

void
PDM_rotation_homogeneous_matrix_to_axis_angle
(
  const double     *homogeneous_matrix,
  const PDM_bool_t  reverse,
        double      axis[3],
        double     *angle
);

/**
 *
 * \brief Converts a rotation expressed as a 4-by-4 homogeneous rotation matrix to euler angles
 *
 * \param [in]   homogeneous_matrix 4-by-4 homogeneous rotation matrix
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [in]   order              Order of rotations to apply
 * \param [in]   intrinsic          <a href="https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations">Axis conventions</a>
 * \param [out]  ang_x              Rotation angle around the x-axis
 * \param [out]  ang_y              Rotation angle around the y-axis
 * \param [out]  ang_z              Rotation angle around the z-axis
 *
 */

void
PDM_rotation_homogeneous_matrix_to_euler_angles
(
  const double     *homogeneous_matrix,
  const PDM_bool_t  reverse,
  const int         order[3],
  const PDM_bool_t  intrinsic,
        double     *ang_x,
        double     *ang_y,
        double     *ang_z
);

/**
 *
 * \brief Converts a rotation expressed as a 4-by-4 homogeneous matrix to a 3-by-3 rotation matrix
 *
 * \param [in]   homogeneous_matrix 4-by-4 rotation matrix
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  rotation_matrix    3-by-3 rotation matrix
 *
 */

void
PDM_rotation_homogeneous_matrix_to_rotation_matrix
(
  const double     *homogeneous_matrix,
  const PDM_bool_t  reverse,
        double     *rotation_matrix
);

// 2 unit vectors to other formats ---

/**
 *
 * \brief Converts a rotation expressed as a rotation from one unit vector to another to an axis-angle representation
 *
 * \param [in]   vector_1 First 3D-vector
 * \param [in]   vector_2 Second 3D-vector
 * \param [in]   reverse  If True, encodes the reverse transformation
 * \param [out]  axis     Rotation axis (3D (unit) vector)
 * \param [out]  angle    Rotation angle (in radians)
 *
 */

void
PDM_rotation_two_vectors_to_axis_angle
(
  const double      vector_1[3],
  const double      vector_2[3],
  const PDM_bool_t  reverse,
        double      axis[3],
        double     *angle
);

/**
 *
 * \brief Converts a rotation expressed as a rotation from one unit vector to another to euler angles
 *
 * \param [in]   vector_1    First 3D-vector
 * \param [in]   vector_2    Second 3D-vector
 * \param [in]   reverse     If True, encodes the reverse transformation
 * \param [in]   order       Order of rotations to apply
 * \param [in]   intrinsic   <a href="https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations">Axis conventions</a>
 * \param [out]  ang_x       Rotation angle around the x-axis
 * \param [out]  ang_y       Rotation angle around the y-axis
 * \param [out]  ang_z       Rotation angle around the z-axis
 *
 */

void
PDM_rotation_two_vectors_to_euler_angles
(
  const double      vector_1[3],
  const double      vector_2[3],
  const PDM_bool_t  reverse,
  const int         order[3],
  const PDM_bool_t  intrinsic,
        double     *ang_x,
        double     *ang_y,
        double     *ang_z
);

/**
 *
 * \brief Converts a rotation expressed as a rotation from one unit vector to another to a 3-by-3 rotation matrix
 *
 * \param [in]   vector_1         First 3D-vector
 * \param [in]   vector_2         Second 3D-vector
 * \param [in]   reverse          If True, encodes the reverse transformation
 * \param [out]  rotation_matrix  3-by-3 rotation matrix
 *
 */

void
PDM_rotation_two_vectors_to_rotation_matrix
(
  const double      vector_1[3],
  const double      vector_2[3],
  const PDM_bool_t  reverse,
        double     *rotation_matrix
);

/**
 *
 * \brief Converts a rotation expressed as a rotation from one unit vector to another to a 4-by-4 homogeneous rotation matrix
 *
 * \param [in]   vector_1           First 3D-vector
 * \param [in]   vector_2           Second 3D-vector
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  homogeneous_matrix 4-by-4 homogeneous rotation matrix
 *
 */

void
PDM_rotation_two_vectors_to_homogeneous_matrix
(
  const double      vector_1[3],
  const double      vector_2[3],
  const PDM_bool_t  reverse,
        double     *homogeneous_matrix
);

/**
 *
 * \brief Converts a rotation expressed as a rotation from one unit vector to another and a rotation center to a 4-by-4 homogeneous rotation matrix
 *
 * \param [in]   vector_1           First 3D-vector
 * \param [in]   vector_2           Second 3D-vector
 * \param [in]   rotation_center    Rotation center
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  homogeneous_matrix 4-by-4 homogeneous rotation matrix
 *
 */

void
PDM_rotation_two_vectors_and_rotation_center_to_homogeneous_matrix
(
  const double      vector_1[3],
  const double      vector_2[3],
  const double      rotation_center[3],
  const PDM_bool_t  reverse,
        double     *homogeneous_matrix
);

// Other formats ---

/**
 *
 * \brief Computes the homogeneous matrix corresponding to switch from cartesian coordinate system A to system B. Axes and origin arguments describe the output coordinate system B using the input coordinate system A.
 *
 * \param [in]   axis_1             First axis 3D-vector
 * \param [in]   axis_2             Second axis 3D-vector
 * \param [in]   axis_3             Thrid axis 3D-vector
 * \param [in]   origin             Origin 3D-vector
 * \param [in]   reverse            If True, encodes the reverse transformation
 * \param [out]  homogeneous_matrix 4-by-4 homogeneous rotation matrix
 *
 */

void
PDM_rotation_axes_and_origin_to_homogeneous_matrix
(
  const double      axis_1[3],
  const double      axis_2[3],
  const double      axis_3[3],
  const double      origin[3],
  const PDM_bool_t  reverse,
        double     *homogeneous_matrix
);

/*----------------------------------------------------------------------------
 *  PDM_rotation_t HOMOGENOUS MATRICES
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Computes A*B -> C
 *
 * \param [in]   A  n-by-n matrix (Row major (C-order) flatten)
 * \param [in]   B  n-by-n matrix (Row major (C-order) flatten)
 * \param [in]   n  size of the square matrices
 * \param [out]  C  n-by-n matrix (Row major (C-order) flatten)
 *
 */


void
PDM_rotation_multiply_n_by_n_matrices
(
  const double *A,
  const double *B,
  const int     n,
        double *C
);

/**
 *
 * \brief Computes A*x -> y
 *
 * \param [in]   A      n-by-n matrix (Row major (C-order) flatten)
 * \param [in]   x      n-by-n_samp vector (Row major (C-order) flatten)
 * \param [out]  y      n-by-n_samp output vector (Row major (C-order) flatten)
 *
 */

void
PDM_rotation_apply_n_by_n_matrix
(
  const double *A,
  const double *x,
  const int     n,
  const int     n_samp,
        double *y
);


/**
 *
 * \brief Applies an homogeneous (4-by-4) matrix to (a) vector(s)
 *
 * \param [in]   homogeneous_matrix 4-by-4 homogeneous matrix (Row major (C-order) flatten)
 * \param [in]   vector             n-by-n_samp vector (Row major (C-order) flatten)
 * \param [in]   n_samp             Number of vector samples
 * \param [out]  vector_out         n-by-n_samp output vector (Row major (C-order) flatten)
 *
 */

void
PDM_rotation_apply_homogeneous_matrix
(
  const double  homogeneous_matrix[16],
  const double *vector,
  const int     n_samp,
        double *vector_out
);

/**
 *
 * \brief Composes (multipies) homogeneous (4-by-4) matrices
 *
 * Computes \f$A_0 \times A_1 \times ... \times A_{n-1}\f$
 *
 * \param [in]   homogeneous_matrices Array of 4-by-4 homogeneous matrix (Row major (C-order) flatten)
 * \param [in]   n_matrices           Number of matrices
 * \param [out]  output_matrix        The resulting 4-by-4 homogeneous matrix
 *
 */

void
PDM_rotation_compose_homogeneous_matrices
(
  const double **homogeneous_matrices,
  const int      n_matrices,
        double   output_matrix[16]
);


/*----------------------------------------------------------------------------
 *  PDM_rotation_t COMPOSITE FUNCTIONS
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Applies the rotation corresponding to the euler angles around the provided rotation center to (a) vector(s)
 *
 * \param [in]   ang_x            Rotation angle around the x-axis
 * \param [in]   ang_y            Rotation angle around the y-axis
 * \param [in]   ang_z            Rotation angle around the z-axis
 * \param [in]   order            Order of rotations to apply
 * \param [in]   intrinsic        <a href="https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations">Axis conventions</a>
 * \param [in]   rotation_center  3D rotation center
 * \param [in]   reverse          If True, encodes the reverse transformation
 * \param [in]   vector           n-by-n_samp vector (Row major (C-order) flatten)
 * \param [in]   n_samp           Number of vector samples
 * \param [out]  vector_out       n-by-n_samp output vector (Row major (C-order) flatten)
 *
 */
void
PDM_rotation_apply_euler_angles_and_rotation_center
(
  const double      ang_x,
  const double      ang_y,
  const double      ang_z,
  const int         order[3],
  const PDM_bool_t  intrinsic,
  const double      rotation_center[3],
  const PDM_bool_t  reverse,
  const double     *vector,
  const int         n_samp,
        double     *vector_out
);

/**
 *
 * \brief Applies the rotation corresponding to the axis and angle around the provided rotation center to (a) vector(s)
 *
 * \param [in]   axis             Rotation axis (3D (unit) vector)
 * \param [in]   angle            Rotation angle (in radians)
 * \param [in]   rotation_center  3D rotation center
 * \param [in]   reverse          If True, encodes the reverse transformation
 * \param [in]   vector           n-by-n_samp vector (Row major (C-order) flatten)
 * \param [in]   n_samp           Number of vector samples
 * \param [out]  vector_out       n-by-n_samp output vector (Row major (C-order) flatten)
 *
 */

void
PDM_rotation_apply_axis_angle_and_rotation_center
(
  const double      axis[3],
  const double      angle,
  const double      rotation_center[3],
  const PDM_bool_t  reverse,
  const double     *vector,
  const int         n_samp,
        double     *vector_out
);

/**
 *
 * \brief Applies the rotation corresponding to the 3-by-3 rotation matrix around the provided rotation center to (a) vector(s)
 *
 * \param [in]   rotation_matrix  3-by-3 rotation matrix
 * \param [in]   rotation_center  3D rotation center
 * \param [in]   reverse          If True, encodes the reverse transformation
 * \param [in]   vector           n-by-n_samp vector (Row major (C-order) flatten)
 * \param [in]   n_samp           Number of vector samples
 * \param [out]  vector_out       n-by-n_samp output vector (Row major (C-order) flatten)
 *
 */

void
PDM_rotation_apply_rotation_matrix_and_rotation_center
(
  const double      rotation_matrix[9],
  const double      rotation_center[3],
  const PDM_bool_t  reverse,
  const double     *vector,
  const int         n_samp,
        double     *vector_out
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_ROTATION_H__ */
