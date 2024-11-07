cdef extern from "pdm_quaternion.h":
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  void PDM_quaternion_axis_angle_to_euler_angles(const double     axis[3],
                                                 const double     angle,
                                                 const int        order[3],
                                                 const PDM_bool_t intrinsic,
                                                       double*    ang_x,
                                                       double*    ang_y,
                                                       double*    ang_z)

  void PDM_quaternion_axis_angle_to_rotation_matrix(const double  axis[3],
                                                    const double  angle,
                                                          double* rotation_matrix)
  void PDM_quaternion_axis_angle_to_homogeneous_matrix(const double  axis[3],
                                                       const double  angle,
                                                             double* homogeneous_matrix)
  void PDM_quaternion_euler_angles_to_axis_angle(const double     ang_x,
                                                 const double     ang_y,
                                                 const double     ang_z,
                                                 const int        order[3],
                                                      PDM_bool_t  intrinsic,
                                                      double      axis[3],
                                                      double*     angle)
  void PDM_quaternion_euler_angles_to_euler_angles(const double     input_ang_x,
                                                   const double     input_ang_y,
                                                   const double     input_ang_z,
                                                   const int        input_order[3],
                                                   const PDM_bool_t input_intrinsic,
                                                   const int        output_order[3],
                                                   const PDM_bool_t output_intrinsic,
                                                         double*    output_ang_x,
                                                         double*    output_ang_y,
                                                         double*    output_ang_z)

  void PDM_quaternion_euler_angles_to_rotation_matrix(const double      ang_x,
                                                      const double      ang_y,
                                                      const double      ang_z,
                                                      const int         order[3],
                                                            PDM_bool_t  intrinsic,
                                                            double*     rotation_matrix)

  void PDM_quaternion_euler_angles_to_homogeneous_matrix(const double     ang_x,
                                                         const double     ang_y,
                                                         const double     ang_z,
                                                         const int        order[3],
                                                               PDM_bool_t intrinsic,
                                                               double*    homogeneous_matrix)

  void PDM_quaternion_rotation_matrix_to_axis_angle(const double* rotation_matrix,
                                                          double  axis[3],
                                                          double* angle)

  void PDM_quaternion_rotation_matrix_to_euler_angles(const double*    rotation_matrix,
                                                      const int        order[3],
                                                      const PDM_bool_t intrinsic,
                                                            double*    ang_x,
                                                            double*    ang_y,
                                                            double*    ang_z)

  void PDM_quaternion_rotation_matrix_to_homogeneous_matrix(const double* rotation_matrix,
                                                                  double* homogeneous_matrix)

  void PDM_quaternion_homogeneous_matrix_to_axis_angle(const double* homogeneous_matrix,
                                                             double  axis[3],
                                                             double* angle)

  void PDM_quaternion_homogeneous_matrix_to_euler_angles(const double*    homogeneous_matrix,
                                                         const int        order[3],
                                                         const PDM_bool_t intrinsic,
                                                               double*    ang_x,
                                                               double*    ang_y,
                                                               double*    ang_z)

  void PDM_quaternion_homogeneous_matrix_to_rotation_matrix(const double* homogeneous_matrix,
                                                                  double* rotation_matrix)

  void PDM_quaternion_two_vectors_to_axis_angle(const double  vector_1[3],
                                                const double  vector_2[3],
                                                      double  axis[3],
                                                      double* angle)

  void PDM_quaternion_two_vectors_to_euler_angles(const double     vector_1[3],
                                                  const double     vector_2[3],
                                                  const int        order[3],
                                                  const PDM_bool_t intrinsic,
                                                        double*    ang_x,
                                                        double*    ang_y,
                                                        double*    ang_z)

  void PDM_quaternion_two_vectors_to_rotation_matrix(const double  vector_1[3],
                                                     const double  vector_2[3],
                                                           double* rotation_matrix)

  void PDM_quaternion_two_vectors_to_homogeneous_matrix(const double  vector_1[3],
                                                        const double  vector_2[3],
                                                              double* homogeneous_matrix)

  void PDM_quaternion_identity_to_homogeneous_matrix(double* homogeneous_matrix)

  void PDM_quaternion_multiply_n_by_n_matrices(const double* A,
                                               const double* B,
                                               const int     n,
                                                     double* C)

  void PDM_quaternion_apply_n_by_n_matrix(const double* A,
                                          const double* x,
                                          const int     n,
                                          const int     n_samp,
                                                double* y_out)

  void PDM_quaternion_apply_translation(const double  translation_vector[3],
                                        const double* vector,
                                        const int     n_samp,
                                              double* vector_out)

  void PDM_quaternion_apply_homogeneous_matrix(const double  homogeneous_matrix[16],
                                               const double* vector,
                                               const int     n_samp,
                                                     double* vector_out)

  void PDM_quaternion_compose_homogeneous_matrices(const double** homogeneous_matrices,
                                                   const int      n_matrices,
                                                         double   output_matrix[16])

  void PDM_quaternion_translation_to_homogeneous_matrix(const double     translation_vector[3],
                                                              PDM_bool_t reverse,
                                                              double     homogeneous_matrix[16])

  void PDM_quaternion_apply_euler_angles_and_rotation_center(const double     ang_x,
                                                             const double     ang_y,
                                                             const double     ang_z,
                                                             const int        order[3],
                                                             const PDM_bool_t intrinsic,
                                                             const double     rotation_center[3],
                                                             const PDM_bool_t reverse,
                                                             const double*    vector,
                                                             const int        n_samp,
                                                                   double*    vector_out)

  void PDM_quaternion_apply_axis_angle_and_rotation_center(const double     axis[3],
                                                           const double     angle,
                                                           const double     rotation_center[3],
                                                           const PDM_bool_t reverse,
                                                           const double*    vector,
                                                           const int        n_samp,
                                                                 double*    vector_out)

  void PDM_quaternion_apply_rotation_matrix_and_rotation_center(const double     rotation_matrix[9],
                                                                const double     rotation_center[3],
                                                                const PDM_bool_t reverse,
                                                                const double*    vector,
                                                                const int        n_samp,
                                                                      double*    vector_out)

# NPY.double_t
# int
# create_numpy_d

def apply_euler_angles_and_rotation_center_to_coords(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] coords,
    NPY.double_t ang_x,
    NPY.double_t ang_y,
    NPY.double_t ang_z,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center = NPY.array([0.,0.,0.],dtype=NPY.double),
    bint reverse = False):
  """
  apply_euler_angles_and_rotation_center_to_coords(coords,ang_x,ang_y,ang_z,order=[2,1,0],intrinsic=True,rotation_center=[0.,0.,0.],reverse=False)

  Applies the rigid transform corresponding to the euler angles around the
  rotation center to the provided coords array

  Parameters:
    coords          (np.ndarray[np.double_t]) : Vector of coordinates (shape = (*n*,3)))
    ang_x           (double)                  : Rotation angle around the x-axis
    ang_y           (double)                  : Rotation angle around the y-axis
    ang_z           (double)                  : Rotation angle around the z-axis
    order           (np.ndarray[np.int32_t])  : Order of rotations to apply 
    intrinsic       (bool)                    : Axis conventions (https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations)
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True, applies the reverse transformation
  """
  # checking sizes
  if order.size != 3:
    raise AssertionError(f"'order' argument of invalid shape {NPY.shape(order)}, expects (3,)")
  if rotation_center.size != 3:
    raise AssertionError(f"'rotation_center' argument of invalid shape {NPY.shape(rotation_center)}, expects (3,)")
  if coords.shape[1] != 3:
    raise AssertionError(f"'coords' argument of invalid shape {NPY.shape(coords)}, expects (n,3)")
  cdef int* order_data = np_to_int_pointer(order)
  cdef double* rotation_center_data = np_to_double_pointer(rotation_center)
  cdef int n_samp = coords.shape[0]
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_out = NPY.empty((n_samp,3),dtype=NPY.double)
  PDM_quaternion_apply_euler_angles_and_rotation_center(ang_x,
                                                        ang_y,
                                                        ang_z,
                                                        order_data,
                                                        <PDM_bool_t> intrinsic,
                                                        rotation_center_data,
                                                        <PDM_bool_t> reverse,
                                                        <double*> coords.data,
                                                        n_samp,
                                                        <double*> vector_out.data)
  return vector_out
  
def apply_euler_angles_and_rotation_center_to_vector_field(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_field,
    NPY.double_t ang_x,
    NPY.double_t ang_y,
    NPY.double_t ang_z,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center = NPY.array([0.,0.,0.],dtype=NPY.double),
    bint reverse = False):
  """
  apply_euler_angles_and_rotation_center_to_vector_field(vector_field,ang_x,ang_y,ang_z,order=[2,1,0],intrinsic=True,rotation_center=[0.,0.,0.],reverse=False)

  Applies the rigid transform corresponding to the euler angles around the
  rotation center to the provided vector field array

  Parameters:
    vector_field    (np.ndarray[np.double_t]) : Vector (shape = (*n*,3)))
    ang_x           (double)                  : Rotation angle around the x-axis
    ang_y           (double)                  : Rotation angle around the y-axis
    ang_z           (double)                  : Rotation angle around the z-axis
    order           (np.ndarray[np.int32_t])  : Order of rotations to apply 
    intrinsic       (bool)                    : Axis conventions (https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations)
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True, applies the reverse transformation
  """
    # checking sizes
  if order.size != 3:
    raise AssertionError(f"'order' argument of invalid shape {NPY.shape(order)}, expects (3,)")
  if rotation_center.size != 3:
    raise AssertionError(f"'rotation_center' argument of invalid shape {NPY.shape(rotation_center)}, expects (3,)")
  if vector_field.shape[1] != 3:
    raise AssertionError(f"'vector_field' argument of invalid shape {NPY.shape(vector_field)}, expects (n,3)")
  cdef int* order_data = np_to_int_pointer(order)
  cdef double rotation_center_data[3]
  rotation_center_data[:] = [0.,0.,0.]
  cdef int n_samp = vector_field.shape[0]
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_out = NPY.empty((n_samp,3),dtype=NPY.double)
  PDM_quaternion_apply_euler_angles_and_rotation_center(ang_x,
                                                        ang_y,
                                                        ang_z,
                                                        order_data,
                                                        <PDM_bool_t> intrinsic,
                                                        rotation_center_data,
                                                        <PDM_bool_t> reverse,
                                                        <double*> vector_field.data,
                                                        n_samp,
                                                        <double*> vector_out.data)
  return vector_out

def apply_axis_angle_and_rotation_center():
  """
  .. todo:: to be implemented
  """
  pass
def apply_rotation_matrix_and_rotation_center():
  """
  .. todo:: to be implemented
  """
  pass