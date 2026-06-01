cdef extern from "pdm_rotation.h":
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # axis angle -> other formats ---

  void PDM_rotation_axis_angle_to_euler_angles(const double     axis[3],
                                               const double     angle,
                                               const PDM_bool_t reverse,
                                               const int        order[3],
                                               const PDM_bool_t intrinsic,
                                                     double*    ang_x,
                                                     double*    ang_y,
                                                     double*    ang_z)

  void PDM_rotation_axis_angle_to_rotation_matrix(const double     axis[3],
                                                  const double     angle,
                                                  const PDM_bool_t reverse,
                                                        double*    rotation_matrix)

  void PDM_rotation_axis_angle_to_homogeneous_matrix(const double     axis[3],
                                                     const double     angle,
                                                     const PDM_bool_t reverse,
                                                           double*    homogeneous_matrix)

  void PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix(const double     axis[3],
                                                                         const double     angle,
                                                                         const double     rotation_center[3],
                                                                         const PDM_bool_t reverse,
                                                                               double*    homogeneous_matrix)

  # Euler angles -> other formats ---

  void PDM_rotation_euler_angles_to_axis_angle(const double     ang_x,
                                               const double     ang_y,
                                               const double     ang_z,
                                               const int        order[3],
                                               const PDM_bool_t intrinsic,
                                               const PDM_bool_t reverse,
                                                     double     axis[3],
                                                     double*    angle)

  void PDM_rotation_euler_angles_to_euler_angles(const double     input_ang_x,
                                                 const double     input_ang_y,
                                                 const double     input_ang_z,
                                                 const int        input_order[3],
                                                 const PDM_bool_t input_intrinsic,
                                                 const PDM_bool_t reverse,
                                                 const int        output_order[3],
                                                 const PDM_bool_t output_intrinsic,
                                                       double*    output_ang_x,
                                                       double*    output_ang_y,
                                                       double*    output_ang_z)

  void PDM_rotation_euler_angles_to_rotation_matrix(const double      ang_x,
                                                    const double      ang_y,
                                                    const double      ang_z,
                                                    const int         order[3],
                                                    const PDM_bool_t  intrinsic,
                                                    const PDM_bool_t reverse,
                                                          double*     rotation_matrix)

  void PDM_rotation_euler_angles_to_homogeneous_matrix(const double     ang_x,
                                                       const double     ang_y,
                                                       const double     ang_z,
                                                       const int        order[3],
                                                       const PDM_bool_t intrinsic,
                                                       const PDM_bool_t reverse,
                                                             double*    homogeneous_matrix)

  void PDM_rotation_euler_angles_and_rotation_center_to_homogeneous_matrix(const double     ang_x,
                                                                           const double     ang_y,
                                                                           const double     ang_z,
                                                                           const int        order[3],
                                                                           const PDM_bool_t intrinsic,
                                                                           const double     rotation_center[3],
                                                                           const PDM_bool_t reverse,
                                                                                 double*    homogeneous_matrix)

  void PDM_rotation_periodic_t_info_to_homogeneous_matrix(const double rotation_center[3],
                                                          const double rotation_angle[3],
                                                          const double translation[3],
                                                          const PDM_bool_t reverse,
                                                                double* homogeneous_matrix)

  # rotation matrix -> other formats ---

  void PDM_rotation_rotation_matrix_to_axis_angle(const double*    rotation_matrix,
                                                  const PDM_bool_t reverse,
                                                        double     axis[3],
                                                        double*    angle)

  void PDM_rotation_rotation_matrix_to_euler_angles(const double*    rotation_matrix,
                                                    const PDM_bool_t reverse,
                                                    const int        order[3],
                                                    const PDM_bool_t intrinsic,
                                                          double*    ang_x,
                                                          double*    ang_y,
                                                          double*    ang_z)

  void PDM_rotation_rotation_matrix_to_homogeneous_matrix(const double*    rotation_matrix,
                                                          const PDM_bool_t reverse,
                                                                double*    homogeneous_matrix)

  void PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix(const double*    rotation_matrix,
                                                                              const double     rotation_center[3],
                                                                              const PDM_bool_t reverse,
                                                                                    double*    homogeneous_matrix)

  # homogeneous matrix -> other formats ---

  void PDM_rotation_homogeneous_matrix_to_axis_angle(const double* homogeneous_matrix,
                                                     const PDM_bool_t reverse,
                                                           double  axis[3],
                                                           double* angle)

  void PDM_rotation_homogeneous_matrix_to_euler_angles(const double*    homogeneous_matrix,
                                                       const PDM_bool_t reverse,
                                                       const int        order[3],
                                                       const PDM_bool_t intrinsic,
                                                             double*    ang_x,
                                                             double*    ang_y,
                                                             double*    ang_z)

  void PDM_rotation_homogeneous_matrix_to_euler_angles_and_translation(const double*    homogeneous_matrix,
                                                                           const PDM_bool_t reverse,
                                                                           const int        order[3],
                                                                           const PDM_bool_t intrinsic,
                                                                                 double*    ang_x,
                                                                                 double*    ang_y,
                                                                                 double*    ang_z,
                                                                                 double     rotation_center[3])

  void PDM_rotation_homogeneous_matrix_to_periodic_t_info(const double*    homogeneous_matrix,
                                                          const PDM_bool_t reverse,
                                                          const PDM_bool_t compute_rotation_center,
                                                                double     rotation_center[3],
                                                                double     rotation_angle[3],
                                                                double     translation[3])

  void PDM_rotation_homogeneous_matrix_to_rotation_matrix(const double*    homogeneous_matrix,
                                                          const PDM_bool_t reverse,
                                                                double*    rotation_matrix)

  # 2 unit vectors -> other formats ---

  void PDM_rotation_two_vectors_to_axis_angle(const double     vector_1[3],
                                              const double     vector_2[3],
                                              const PDM_bool_t reverse,
                                                    double     axis[3],
                                                    double*    angle)

  void PDM_rotation_two_vectors_to_euler_angles(const double     vector_1[3],
                                                const double     vector_2[3],
                                                const PDM_bool_t reverse,
                                                const int        order[3],
                                                const PDM_bool_t intrinsic,
                                                      double*    ang_x,
                                                      double*    ang_y,
                                                      double*    ang_z)

  void PDM_rotation_two_vectors_to_rotation_matrix(const double     vector_1[3],
                                                   const double     vector_2[3],
                                                   const PDM_bool_t reverse,
                                                         double*    rotation_matrix)

  void PDM_rotation_two_vectors_to_homogeneous_matrix(const double     vector_1[3],
                                                      const double     vector_2[3],
                                                      const PDM_bool_t reverse,
                                                            double*    homogeneous_matrix)

  void PDM_rotation_two_vectors_and_rotation_center_to_homogeneous_matrix(const double     vector_1[3],
                                                                          const double     vector_2[3],
                                                                          const double     rotation_center[3],
                                                                          const PDM_bool_t reverse,
                                                                                double*    homogeneous_matrix)

  # other formats ---

  void PDM_rotation_axes_and_origin_to_homogeneous_matrix(const double     axis_1[3],
                                                          const double     axis_2[3],
                                                          const double     axis_3[3],
                                                          const double     origin[3],
                                                          const PDM_bool_t reverse,
                                                                double*    homogeneous_matrix)
  # (homogeneous) matrix manipulation ---

  void PDM_rotation_apply_homogeneous_matrix(const double  homogeneous_matrix[16],
                                             const double* vector,
                                             const int     n_samp,
                                                   double* vector_out)

  void PDM_rotation_compose_homogeneous_matrices(const double** homogeneous_matrices,
                                                 const int      n_matrices,
                                                       double   output_matrix[16])

  # apply functions ---

  void PDM_rotation_apply_euler_angles_and_rotation_center(const double     ang_x,
                                                           const double     ang_y,
                                                           const double     ang_z,
                                                           const int        order[3],
                                                           const PDM_bool_t intrinsic,
                                                           const double     rotation_center[3],
                                                           const PDM_bool_t reverse,
                                                           const double*    vector,
                                                           const int        n_samp,
                                                                 double*    vector_out)

  void PDM_rotation_apply_axis_angle_and_rotation_center(const double     axis[3],
                                                         const double     angle,
                                                         const double     rotation_center[3],
                                                         const PDM_bool_t reverse,
                                                         const double*    vector,
                                                         const int        n_samp,
                                                               double*    vector_out)

  void PDM_rotation_apply_rotation_matrix_and_rotation_center(const double     rotation_matrix[9],
                                                              const double     rotation_center[3],
                                                              const PDM_bool_t reverse,
                                                              const double*    vector,
                                                              const int        n_samp,
                                                                    double*    vector_out)

# default values ---
# 'buffer' defaults are defined here and then referenced in functions signatures
# to avoid memory leak in sanitize mode (cf MR!100)
_default_order_ = NPY.array([2,1,0],dtype=NPY.int32)
_default_rotation_center_ = NPY.array([0.,0.,0.],dtype=NPY.double)
# ---
# cdef NPY.ndarray[NPY.int32_t, mode='c', ndim=1] _default_order_ = NPY.array([2,1,0],dtype=NPY.int32)

def _check_euler_angles_order(
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order,
    str name = "order"):
  if order.size != 3:
    raise AssertionError(f"'{name}' argument of invalid shape {NPY.shape(order)}, expects (3,)")
  if NPY.any(NPY.sort(order) != NPY.array([0,1,2],dtype=NPY.int32)):
    raise AssertionError(f"'{name}' argument ({order}), expects a permutation of (0,1,2)")

def _check_matrix_shape(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] matrix,
    tuple expec_shape,
    str name = "rotation_matrix"):
  if NPY.shape(matrix) != expec_shape:
    raise AssertionError(f"'{name}' of invalid shape {NPY.shape(matrix)}. Expects {expec_shape}.")

def _check_axis_shape(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis,
    str name = "axis"):
  if axis.size != 3:
    raise AssertionError(f"'{name}' argument of invalid shape {NPY.shape(axis)}, expects (3,)")

# region Axis angle to other formats -------------------------------------------

def axis_angle_to_euler_angles(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis,
    NPY.double_t angle,
    bint reverse = False,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True):
  """axis_angle_to_euler_angles(axis,angle,reverse=False,order=(2,1,0),intrinsic=True)

  Converts a rotation expressed as axis-angle to Euler angles

  Parameters:
    axis            (np.ndarray[np.double_t]) : Rotation axis (shape = (3,))
    angles          (double)                  : Rotation angle (in *radians*)
    reverse         (bool)                    : If True computes the reverse transformation
    order           (np.ndarray[np.int32_t])  : Order of rotations to apply
    intrinsic       (bool)                    : `Axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_

  Returns:
    - Rotation angle around the x-axis (`double`)
    - Rotation angle around the y-axis (`double`)
    - Rotation angle around the z-axis (`double`)
  """
  _check_axis_shape(axis)
  cdef NPY.double_t ang_x,ang_y,ang_z
  cdef int* order_data = np_to_int_pointer(order)
  PDM_rotation_axis_angle_to_euler_angles(<double*>axis.data,
                                          angle,
                                          <PDM_bool_t> reverse,
                                          order_data,
                                          <PDM_bool_t> intrinsic,
                                          &ang_x,
                                          &ang_y,
                                          &ang_z)
  return ang_x,ang_y,ang_z

def axis_angle_to_rotation_matrix(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis,
    NPY.double_t angle,
    bint reverse = False):
  """axis_angle_to_rotation_matrix(axis,angle,reverse=False)

  Converts a rotation expressed as axis-angle to a 3-by-3 rotation matrix

  Parameters:
    axis            (np.ndarray[np.double_t]) : Rotation axis (shape = (3,))
    angles          (double)                  : Rotation angle (in *radians*)
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    3-by-3 rotation matrix (`np.ndarray[np.double_t]`, shape = (3,3))
  """
  _check_axis_shape(axis)
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] rotation_matrix = NPY.empty((3,3),dtype=NPY.double)
  PDM_rotation_axis_angle_to_rotation_matrix(<double*>axis.data,
                                             angle,
                                             <PDM_bool_t> reverse,
                                             <double*>rotation_matrix.data)
  return rotation_matrix

def axis_angle_to_homogeneous_matrix(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis,
    NPY.double_t angle,
    bint reverse = False):
  """axis_angle_to_homogeneous_matrix(axis,angle,reverse=False)

  Converts a rotation expressed as axis-angle to a 4-by-4 homogeneous matrix

  Parameters:
    axis            (np.ndarray[np.double_t]) : Rotation axis (shape = (3,))
    angles          (double)                  : Rotation angle (in *radians*)
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    4-by-4 homogeneous matrix (`np.ndarray[np.double_t]`, shape = (4,4))
  """
  _check_axis_shape(axis)
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix = NPY.empty((4,4),dtype=NPY.double)
  PDM_rotation_axis_angle_to_homogeneous_matrix(<double*>axis.data,
                                                angle,
                                                <PDM_bool_t> reverse,
                                                <double*>homogeneous_matrix.data)
  return homogeneous_matrix

def axis_angle_and_rotation_center_to_homogeneous_matrix(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis,
    NPY.double_t angle,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center=_default_rotation_center_,
    bint reverse = False):
  """axis_angle_and_rotation_center_to_homogeneous_matrix(axis,angle,rotation_center=[0.,0.,0.],reverse=False)

  Computes the homogeneous matrix corresponding to the rotation of the provided angle around the
  provided axis and rotation center

  Parameters:
    axis            (np.ndarray[np.double_t]) : Rotation axis (shape = (3,))
    angles          (double)                  : Rotation angle (in *radians*)
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    4-by-4 homogeneous matrix (`np.ndarray[np.double_t]`, shape = (4,4))
  """
  _check_axis_shape(axis)
  _check_axis_shape(rotation_center,"rotation_center")
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix = NPY.empty((4,4),dtype=NPY.double)
  PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix(<double*>axis.data,
                                                                    angle,
                                                                    <double*>rotation_center.data,
                                                                    <PDM_bool_t> reverse,
                                                                    <double*>homogeneous_matrix.data)
  return homogeneous_matrix

# region Euler angles to other formats -----------------------------------------

def euler_angles_to_axis_angle(
    NPY.double_t ang_x,
    NPY.double_t ang_y,
    NPY.double_t ang_z,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True,
    bint reverse = False):
  """euler_angles_to_axis_angle(ang_x,ang_y,ang_z,order=(2,1,0),intrinsic=True,reverse=False)

  Converts a rotation expressed as Euler angles to axis-angle

  Parameters:
    ang_x           (double)                  : Rotation angle around the x-axis
    ang_y           (double)                  : Rotation angle around the y-axis
    ang_z           (double)                  : Rotation angle around the z-axis
    order           (np.ndarray[np.int32_t])  : Order of rotations to apply
    intrinsic       (bool)                    : `Axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    - Rotation axis  (`np.ndarray[np.double_t]`, shape = (3,))
    - Rotation angle (`double`, in *radians*)
  """
  _check_euler_angles_order(order)
  cdef int* order_data = np_to_int_pointer(order)
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis = NPY.empty((3,),dtype=NPY.double)
  cdef NPY.double_t angle
  PDM_rotation_euler_angles_to_axis_angle(ang_x,
                                          ang_y,
                                          ang_z,
                                          order_data,
                                          <PDM_bool_t> intrinsic,
                                          <PDM_bool_t> reverse,
                                          <double*> axis.data,
                                          &angle)
  return axis,angle

def euler_angles_to_euler_angles(
    NPY.double_t ang_x,
    NPY.double_t ang_y,
    NPY.double_t ang_z,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] input_order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint input_intrinsic = True,
    bint reverse = False,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] output_order = _default_order_,# = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint output_intrinsic = True):
  """euler_angles_to_euler_angles(ang_x,ang_y,ang_z,input_order=(2,1,0),input_intrinsic=True,reverse=False,output_order=(2,1,0),output_intrinsic=True)

  Converts a rotation expressed as Euler angles to another Euler angles expression

  Parameters:
    ang_x            (double)                  : Rotation angle around the x-axis
    ang_y            (double)                  : Rotation angle around the y-axis
    ang_z            (double)                  : Rotation angle around the z-axis
    input_order      (np.ndarray[np.int32_t])  : Input order of rotations to apply
    input_intrinsic  (bool)                    : Input `axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_
    reverse          (bool)                    : If True computes the reverse transformation
    output_order     (np.ndarray[np.int32_t])  : Output order of rotations to apply
    output_intrinsic (bool)                    : Output `axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_

  Returns:
    - Rotation angle around the x-axis (`double`)
    - Rotation angle around the y-axis (`double`)
    - Rotation angle around the z-axis (`double`)
  """
  _check_euler_angles_order(input_order, "input_order")
  _check_euler_angles_order(output_order,"output_order")

  cdef NPY.double_t output_ang_x,output_ang_y,output_ang_z
  cdef int* inp_order_data = np_to_int_pointer(input_order)
  cdef int* out_order_data = np_to_int_pointer(output_order)
  PDM_rotation_euler_angles_to_euler_angles(ang_x,
                                            ang_y,
                                            ang_z,
                                            inp_order_data,
                                            <PDM_bool_t> input_intrinsic,
                                            <PDM_bool_t> reverse,
                                            out_order_data,
                                            <PDM_bool_t> output_intrinsic,
                                            &output_ang_x,
                                            &output_ang_y,
                                            &output_ang_z)
  return output_ang_x,output_ang_y,output_ang_z

def euler_angles_to_rotation_matrix(
    NPY.double_t ang_x,
    NPY.double_t ang_y,
    NPY.double_t ang_z,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True,
    bint reverse = False):
  """
  euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z,order=(2,1,0),intrinsic=True,reverse=False)

  Computes the rotation matrix corresponding to the provided Euler angles

  Parameters:
    ang_x           (double)                  : Rotation angle around the x-axis
    ang_y           (double)                  : Rotation angle around the y-axis
    ang_z           (double)                  : Rotation angle around the z-axis
    order           (np.ndarray[np.int32_t])  : Order of rotations to apply
    intrinsic       (bool)                    : `Axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    3-by-3 rotation matrix (`np.ndarray[np.double_t]`, shape = (3,3))
  """

  _check_euler_angles_order(order)
  cdef int* order_data = np_to_int_pointer(order)
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] rotation_matrix = NPY.empty((3,3),dtype=NPY.double)
  PDM_rotation_euler_angles_to_rotation_matrix(ang_x,
                                               ang_y,
                                               ang_z,
                                               order_data,
                                               <PDM_bool_t> intrinsic,
                                               <PDM_bool_t> reverse,
                                               <double*> rotation_matrix.data)
  return rotation_matrix

def euler_angles_to_homogeneous_matrix(
    NPY.double_t ang_x,
    NPY.double_t ang_y,
    NPY.double_t ang_z,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True,
    bint reverse = False):
  """
  euler_angles_to_homogeneous_matrix(ang_x, ang_y, ang_z, order=(2,1,0), intrinsic=True, reverse=False)

  Computes the homogeneous matrix corresponding to the provided Euler angles

  Parameters:
    ang_x           (double)                  : Rotation angle around the x-axis
    ang_y           (double)                  : Rotation angle around the y-axis
    ang_z           (double)                  : Rotation angle around the z-axis
    order           (np.ndarray[np.int32_t])  : Order of rotations to apply
    intrinsic       (bool)                    : `Axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    4-by-4 homogeneous matrix (`np.ndarray[np.double_t]`, shape = (4,4))
  """

  _check_euler_angles_order(order)
  cdef int* order_data = np_to_int_pointer(order)
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix = NPY.empty((4,4),dtype=NPY.double)
  PDM_rotation_euler_angles_to_homogeneous_matrix(ang_x,
                                                  ang_y,
                                                  ang_z,
                                                  order_data,
                                                  <PDM_bool_t> intrinsic,
                                                  <PDM_bool_t> reverse,
                                                  <double*> homogeneous_matrix.data)
  return homogeneous_matrix

def euler_angles_and_rotation_center_to_homogeneous_matrix(
    NPY.double_t ang_x,
    NPY.double_t ang_y,
    NPY.double_t ang_z,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center = _default_rotation_center_,# = NPY.array([0,0,0],dtype=NPY.float64),
    bint reverse = False):
  """
  euler_angles_and_rotation_center_to_homogeneous_matrix(ang_x,ang_y,ang_z,order=(2,1,0),intrinsic=True,rotation_center=[0.,0.,0.],reverse=False)

  Computes the homogeneous matrix corresponding to the provided Euler angles
  and rotation center.

  Parameters:
    ang_x           (double)                  : Rotation angle around the x-axis
    ang_y           (double)                  : Rotation angle around the y-axis
    ang_z           (double)                  : Rotation angle around the z-axis
    order           (np.ndarray[np.int32_t])  : Order of rotations to apply
    intrinsic       (bool)                    : `Axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    4-by-4 homogeneous matrix (`np.ndarray[np.double_t]`, shape = (4,4))
  """

  _check_euler_angles_order(order)
  _check_axis_shape(rotation_center,"rotation_center")
  cdef int* order_data = np_to_int_pointer(order)
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix = NPY.empty((4,4),dtype=NPY.double)
  PDM_rotation_euler_angles_and_rotation_center_to_homogeneous_matrix(ang_x,
                                                                      ang_y,
                                                                      ang_z,
                                                                      order_data,
                                                                      <PDM_bool_t> intrinsic,
                                                                      <double*> rotation_center.data,
                                                                      <PDM_bool_t> reverse,
                                                                      <double*> homogeneous_matrix.data)
  return homogeneous_matrix

def periodic_t_info_to_homogeneous_matrix(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_angle,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] translation,
    bint reverse = False):
  """periodic_t_info_to_homogeneous_matrix(rotation_center,rotation_angle,translation,reverse=False)

  Converts the info of a CGNS Periodic_t node to a 4-by-4 homogeneous matrix.
  Rotation angles are applied as **intrinsic Euler angles** applied in a *(2,1,0)* order.
  Translation is applied after the rotation.

  Parameters:
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    rotation_angle  (np.ndarray[np.double_t]) : Rotation angles around the x, y and z axes (shape = (3,))
    translation     (np.ndarray[np.double_t]) : Translation vector (shape = (3,))
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    4-by-4 homogeneous matrix (`np.ndarray[np.double_t]`, shape = (4,4))
  """
  _check_axis_shape(rotation_center,"rotation_center")
  _check_axis_shape(rotation_angle,"rotation_angle")
  _check_axis_shape(translation,"translation")
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix = NPY.empty((4,4),dtype=NPY.double)
  PDM_rotation_periodic_t_info_to_homogeneous_matrix(<double*> rotation_center.data,
                                                     <double*> rotation_angle.data,
                                                     <double*> translation.data,
                                                     <PDM_bool_t> reverse,
                                                     <double*> homogeneous_matrix.data)
  return homogeneous_matrix

# region rotation matrix to other formats --------------------------------------

def rotation_matrix_to_axis_angle(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] rotation_matrix,
    bint reverse = False):
  """rotation_matrix_to_axis_angle(rotation_matrix,reverse=False)

  Converts a rotation expressed as a rotation matrix to axis-angle

  Parameters:
    rotation_matrix (np.ndarray[np.double_t]) : Rotation matrix (shape = (3,3))
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    Rotation axis  (`np.ndarray[np.double_t]`, shape = (3,))
    Rotation angle (`double`, in *radians*)
  """
  _check_matrix_shape(rotation_matrix,(3,3))
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis = NPY.empty((3,),dtype=NPY.double)
  cdef NPY.double_t angle
  PDM_rotation_rotation_matrix_to_axis_angle(<double*> rotation_matrix.data,
                                             <PDM_bool_t> reverse,
                                             <double*> axis.data,
                                             &angle)
  return axis,angle

def rotation_matrix_to_euler_angles(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] rotation_matrix,
    bint reverse = False,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True):
  """rotation_matrix_to_euler_angles(rotation_matrix,reverse=False,order=(2,1,0),intrinsic=True)

  Converts a rotation expressed as a rotation matrix to Euler angles

  Parameters:
    rotation_matrix (np.ndarray[np.double_t]) : Rotation matrix (shape = (3,3))
    reverse         (bool)                    : If True computes the reverse transformation
    order           (np.ndarray[np.int32_t])  : Order of rotations to apply
    intrinsic       (bool)                    : `Axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_

  Returns:
    Rotation angle around the x-axis (`double`)
    Rotation angle around the y-axis (`double`)
    Rotation angle around the z-axis (`double`)
  """
  _check_matrix_shape(rotation_matrix,(3,3))
  _check_euler_angles_order(order)

  cdef NPY.double_t ang_x,ang_y,ang_z
  cdef int* out_order_data = np_to_int_pointer(order)
  PDM_rotation_rotation_matrix_to_euler_angles(<double*> rotation_matrix.data,
                                               <PDM_bool_t> reverse,
                                               out_order_data,
                                               <PDM_bool_t> intrinsic,
                                               &ang_x,
                                               &ang_y,
                                               &ang_z)
  return ang_x,ang_y,ang_z

def rotation_matrix_to_homogeneous_matrix(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] rotation_matrix,
    bint reverse = False):
  """rotation_matrix_to_homogeneous_matrix(rotation_matrix,reverse=False)

  Converts a rotation expressed as a rotation matrix to an homogeneous matrix

  Parameters:
    rotation_matrix (np.ndarray[np.double_t]) : Rotation matrix (shape = (3,3))
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    4-by-4 homogeneous matrix (`np.ndarray[np.double_t]`, shape = (4,4))
  """
  _check_matrix_shape(rotation_matrix,(3,3))
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix = NPY.empty((4,4),dtype=NPY.double)
  PDM_rotation_rotation_matrix_to_homogeneous_matrix(<double*> rotation_matrix.data,
                                                     <PDM_bool_t> reverse,
                                                     <double*> homogeneous_matrix.data)
  return homogeneous_matrix

def rotation_matrix_and_rotation_center_to_homogeneous_matrix(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] rotation_matrix,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center=_default_rotation_center_,
    bint reverse = False):
  """rotation_matrix_and_rotation_center_to_homogeneous_matrix(rotation_matrix,rotation_center=[0.,0.,0.],reverse=False)

  Computes the homogeneous matrix corresponding to the provided rotation matrix
  and rotation center

  Parameters:
    rotation_matrix (np.ndarray[np.double_t]) : Rotation matrix (shape = (3,3))
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    4-by-4 homogeneous matrix (`np.ndarray[np.double_t]`, shape = (4,4))

  """
  _check_matrix_shape(rotation_matrix,(3,3))
  _check_axis_shape(rotation_center,"rotation_center")
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix = NPY.empty((4,4),dtype=NPY.double)
  PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix(<double*> rotation_matrix.data,
                                                                         <double*> rotation_center.data,
                                                                         <PDM_bool_t> reverse,
                                                                         <double*> homogeneous_matrix.data)
  return homogeneous_matrix

# region homogeneous matrix to other formats -----------------------------------

def homogeneous_matrix_to_axis_angle(
  NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix,
  bint reverse = False,):
  """homogeneous_matrix_to_axis_angle(homogeneous_matrix,reverse=False)

  Converts a rotation expressed as an homogeneous matrix to axis-angle

  Caution: does not consider the affine part of the transformation (rotation center and translation)

  Parameters:
    homogeneous_matrix (np.ndarray[np.double_t]) : Homogeneous matrix (shape = (4,4))
    reverse            (bool)                    : If True computes the reverse transformation

  Returns:
    Rotation axis  (`np.ndarray[np.double_t]`, shape = (3,))
    Rotation angle (`double`, in *radians*)
  """
  _check_matrix_shape(homogeneous_matrix,(4,4))
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis = NPY.empty((3,),dtype=NPY.double)
  cdef NPY.double_t angle
  PDM_rotation_homogeneous_matrix_to_axis_angle(<double*> homogeneous_matrix.data,
                                                <PDM_bool_t> reverse,
                                                <double*> axis.data,
                                                &angle)
  return axis,angle

def homogeneous_matrix_to_euler_angles(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix,
    bint reverse = False,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True):
  """homogeneous_matrix_to_euler_angles(homogeneous_matrix,reverse=False,order=(2,1,0),intrinsic=True)

  Computes a rotation expressed as an homogeneous matrix to Euler angles

  Caution: does not consider the affine part of the transformation (rotation center and translation)

  Parameters:
    homogeneous_matrix (np.ndarray[np.double_t]) : Homogeneous matrix (shape = (4,4))
    reverse            (bool)                    : If True computes the reverse transformation
    order              (np.ndarray[np.int32_t])  : Order of rotations to apply
    intrinsic          (bool)                    : `Axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_

  Returns:
    Rotation angle around the x-axis (`double`)
    Rotation angle around the y-axis (`double`)
    Rotation angle around the z-axis (`double`)
  """
  _check_matrix_shape(homogeneous_matrix,(4,4))
  _check_euler_angles_order(order)
  cdef NPY.double_t ang_x,ang_y,ang_z
  cdef int* out_order_data = np_to_int_pointer(order)
  PDM_rotation_homogeneous_matrix_to_euler_angles(<double*> homogeneous_matrix.data,
                                                  <PDM_bool_t> reverse,
                                                  out_order_data,
                                                  <PDM_bool_t> intrinsic,
                                                  &ang_x,
                                                  &ang_y,
                                                  &ang_z)
  return ang_x,ang_y,ang_z

def homogeneous_matrix_to_euler_angles_and_rotation_center(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix,
    bint reverse = False,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True):
  """homogeneous_matrix_to_euler_angles_and_rotation_center(homogeneous_matrix,reverse=False,order=(2,1,0),intrinsic=True)

  Computes a rotation expressed as an homogeneous matrix to Euler angles and a rotation center

  Caution: does not consider the rotation-axis-wise translation component

  Parameters:
    homogeneous_matrix (np.ndarray[np.double_t]) : Homogeneous matrix (shape = (4,4))
    reverse            (bool)                    : If True computes the reverse transformation
    order              (np.ndarray[np.int32_t])  : Order of rotations to apply
    intrinsic          (bool)                    : `Axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_

  Returns:
    Rotation angle around the x-axis (`double`)
    Rotation angle around the y-axis (`double`)
    Rotation angle around the z-axis (`double`)
    Rotation center (`np.ndarray[np.double_t]`, shape = (3,))
  """
  _check_matrix_shape(homogeneous_matrix,(4,4))
  _check_euler_angles_order(order)
  cdef NPY.double_t ang_x,ang_y,ang_z
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center = NPY.empty((3,),dtype=NPY.double)
  cdef int* out_order_data = np_to_int_pointer(order)
  PDM_rotation_homogeneous_matrix_to_euler_angles_and_translation(<double*> homogeneous_matrix.data,
                                                                      <PDM_bool_t> reverse,
                                                                      out_order_data,
                                                                      <PDM_bool_t> intrinsic,
                                                                      &ang_x,
                                                                      &ang_y,
                                                                      &ang_z,
                                                                      <double*> rotation_center.data)
  return ang_x,ang_y,ang_z,rotation_center

def homogeneous_matrix_to_periodic_t_info(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix,
    bint reverse = False,
    bint compute_rotation_center = False,
    ):
  """homogeneous_matrix_to_periodic_t_info(homogeneous_matrix,reverse=False)

  Converts a 4-by-4 homogeneous matrix to the info of a CGNS Periodic_t node
  Rotation angles are applied as **intrinsic Euler angles** applied in a *(2,1,0)* order.
  Translation is applied after the rotation.

  Parameters:
    homogeneous_matrix      (np.ndarray[np.double_t]) : Homogeneous matrix (shape = (4,4))
    reverse                 (bool)                    : If True computes the reverse transformation
    compute_rotation_center (bool)                    : If True computes the rotation center so that the translation is along the rotation axis, set to 0. otherwise

  Returns:
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    rotation_angle  (np.ndarray[np.double_t]) : Rotation angles around the x, y and z axes (shape = (3,))
    translation     (np.ndarray[np.double_t]) : Translation vector (shape = (3,))
  """
  _check_matrix_shape(homogeneous_matrix,(4,4))
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center = NPY.empty((3,),dtype=NPY.double)
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_angle = NPY.empty((3,),dtype=NPY.double)
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=1] translation = NPY.empty((3,),dtype=NPY.double)
  PDM_rotation_homogeneous_matrix_to_periodic_t_info(<double*> homogeneous_matrix.data,
                                                     <PDM_bool_t> reverse,
                                                     <PDM_bool_t> compute_rotation_center,
                                                     <double*> rotation_center.data,
                                                     <double*> rotation_angle.data,
                                                     <double*> translation.data)
  return rotation_center,rotation_angle,translation

def homogeneous_matrix_to_rotation_matrix(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix,
    bint reverse = False):
  """homogeneous_matrix_to_rotation_matrix(homogeneous_matrix,reverse=False)

  Converts a rotation expressed as an homogeneous matrix to a rotation matrix

  Parameters:
    homogeneous_matrix      (np.ndarray[np.double_t]) : Homogeneous matrix (shape = (4,4))
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    3-by-3 rotation matrix (`np.ndarray[np.double_t]`, shape = (3,3))
  """
  _check_matrix_shape(homogeneous_matrix,(4,4))
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] rotation_matrix = NPY.empty((3,3),dtype=NPY.double)
  PDM_rotation_homogeneous_matrix_to_rotation_matrix(<double*> homogeneous_matrix.data,
                                                     <PDM_bool_t> reverse,
                                                     <double*> rotation_matrix.data)
  return rotation_matrix

# region 2 unit vectors to other formats ---------------------------------------

def two_vectors_to_axis_angle(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] vector_1,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] vector_2,
    bint reverse = False):
  """two_vectors_to_axis_angle(vector_1,vector_2,reverse=False)

  Converts a rotation from the first vector to the latter to a rotation expressed as axis angle

  Parameters:
    vector_1        (np.ndarray[np.double_t]) : First vector  (shape = (3,))
    vector_2        (np.ndarray[np.double_t]) : Second vector (shape = (3,))
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    Rotation axis  (`np.ndarray[np.double_t]`, shape = (3,))
    Rotation angle (`double`, in *radians*)
  """
  _check_axis_shape(vector_1,"vector_1")
  _check_axis_shape(vector_2,"vector_2")
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis = NPY.empty((3,),dtype=NPY.double)
  cdef NPY.double_t angle
  PDM_rotation_two_vectors_to_axis_angle(<double*> vector_1.data,
                                         <double*> vector_2.data,
                                         <PDM_bool_t> reverse,
                                         <double*> axis.data,
                                         &angle)
  return axis,angle

def two_vectors_to_euler_angles(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] vector_1,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] vector_2,
    bint reverse = False,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True):
  """two_vectors_to_euler_angles(vector_1,vector_2,reverse=False,order=(2,1,0),intrinsic=True)

  Converts a rotation from the first vector to the latter to a rotation expressed as Euler angles

  Parameters:
    vector_1        (np.ndarray[np.double_t]) : First vector  (shape = (3,))
    vector_2        (np.ndarray[np.double_t]) : Second vector (shape = (3,))
    reverse         (bool)                    : If True computes the reverse transformation
    order           (np.ndarray[np.int32_t])  : Order of rotations to apply
    intrinsic       (bool)                    : `Axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_

  Returns:
    Rotation angle around the x-axis (`double`)
    Rotation angle around the y-axis (`double`)
    Rotation angle around the z-axis (`double`)

  """
  _check_axis_shape(vector_1,"vector_1")
  _check_axis_shape(vector_2,"vector_2")
  _check_euler_angles_order(order,"order")
  cdef NPY.double_t ang_x,ang_y,ang_z
  cdef int* order_data = np_to_int_pointer(order)
  PDM_rotation_two_vectors_to_euler_angles(<double*> vector_1.data,
                                           <double*> vector_2.data,
                                           <PDM_bool_t> reverse,
                                           order_data,
                                           <PDM_bool_t> intrinsic,
                                           &ang_x,
                                           &ang_y,
                                           &ang_z)
  return ang_x,ang_y,ang_z

def two_vectors_to_rotation_matrix(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] vector_1,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] vector_2,
    bint reverse = False):
  """two_vectors_to_rotation_matrix(vector_1,vector_2,reverse=False)

  Converts a rotation from the first vector to the latter to a rotation expressed as a 3-by-3 rotation matrix

  Parameters:
    vector_1        (np.ndarray[np.double_t]) : First vector  (shape = (3,))
    vector_2        (np.ndarray[np.double_t]) : Second vector (shape = (3,))
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    3-by-3 rotation matrix (`np.ndarray[np.double_t]`, shape = (3,3))
  """
  _check_axis_shape(vector_1,"vector_1")
  _check_axis_shape(vector_2,"vector_2")
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] rotation_matrix = NPY.empty((3,3),dtype=NPY.double)
  PDM_rotation_two_vectors_to_rotation_matrix(<double*> vector_1.data,
                                              <double*> vector_2.data,
                                              <PDM_bool_t> reverse,
                                              <double*> rotation_matrix.data)
  return rotation_matrix

def two_vectors_to_homogeneous_matrix(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] vector_1,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] vector_2,
    bint reverse = False):
  """two_vectors_to_homogeneous_matrix(vector_1,vector_2,reverse=False)

  Converts a rotation from the first vector to the latter to a rotation expressed as a 4-by-4 homogeneous matrix

  Parameters:
    vector_1        (np.ndarray[np.double_t]) : First vector  (shape = (3,))
    vector_2        (np.ndarray[np.double_t]) : Second vector (shape = (3,))
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    4-by-4 homogeneous matrix (`np.ndarray[np.double_t]`, shape = (4,4))
  """
  _check_axis_shape(vector_1,"vector_1")
  _check_axis_shape(vector_2,"vector_2")
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix = NPY.empty((4,4),dtype=NPY.double)
  PDM_rotation_two_vectors_to_homogeneous_matrix(<double*> vector_1.data,
                                                 <double*> vector_2.data,
                                                 <PDM_bool_t> reverse,
                                                 <double*> homogeneous_matrix.data)
  return homogeneous_matrix

def two_vectors_and_rotation_center_to_homogeneous_matrix(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] vector_1,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] vector_2,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center=_default_rotation_center_,
    bint reverse = False):
  """two_vectors_and_rotation_center_to_homogeneous_matrix(vector_1,vector_2,rotation_center=[0.,0.,0.],reverse=False)

  Computes the homogeneous matrix corresponding to a rotation from the first vector to the latter
  around the provided rotation center.

  Parameters:
    vector_1        (np.ndarray[np.double_t]) : First vector  (shape = (3,))
    vector_2        (np.ndarray[np.double_t]) : Second vector (shape = (3,))
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True computes the reverse transformation

  Returns:
    4-by-4 homogeneous matrix (`np.ndarray[np.double_t]`, shape = (4,4))
  """
  _check_axis_shape(vector_1,"vector_1")
  _check_axis_shape(vector_2,"vector_2")
  _check_axis_shape(rotation_center,"rotation_center")
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix = NPY.empty((4,4),dtype=NPY.double)
  PDM_rotation_two_vectors_and_rotation_center_to_homogeneous_matrix(<double*> vector_1.data,
                                                                     <double*> vector_2.data,
                                                                     <double*> rotation_center.data,
                                                                     <PDM_bool_t> reverse,
                                                                     <double*> homogeneous_matrix.data)
  return homogeneous_matrix


# region change in reference frame ---------------------------------------------

def axes_and_origin_to_homogeneous_matrix(
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis_1,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis_2,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis_3,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] origin,
    bint reverse = False):
  """axes_and_origin_to_homogeneous_matrix(axis_1, axis_2, axis_3, origin, reverse=False)

  Computes the homogeneous matrix corresponding to switch from cartesian coordinate
  system A to system B. Axes and origin arguments describe the output coordinate
  system B using the input coordinate system A

  Parameters:
    axis_1          (np.ndarray[np.double_t]) : First axis (shape = (3,))
    axis_2          (np.ndarray[np.double_t]) : Second axis (shape = (3,))
    axis_3          (np.ndarray[np.double_t]) : Third axis (shape = (3,))
    origin          (np.ndarray[np.double_t]) : Origin (shape = (3,))
    reverse         (bool)                    : If True returns the homogeneous matrix corresponding to the inverse transformation

  Returns:
    4-by-4 homogeneous matrix (`np.ndarray[np.double_t]`, shape = (4,4))
  """
  _check_axis_shape(axis_1,"axis_1")
  _check_axis_shape(axis_2,"axis_2")
  _check_axis_shape(axis_3,"axis_3")
  _check_axis_shape(origin,"origin")
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix = NPY.empty((4,4),dtype=NPY.double)
  PDM_rotation_axes_and_origin_to_homogeneous_matrix(<double*> axis_1.data,
                                                     <double*> axis_2.data,
                                                     <double*> axis_3.data,
                                                     <double*> origin.data,
                                                     <PDM_bool_t> reverse,
                                                     <double*> homogeneous_matrix.data)
  return homogeneous_matrix

# region 'Apply' functions -----------------------------------------------------

def apply_euler_angles_and_rotation_center_to_coords(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] coords,
    NPY.double_t ang_x,
    NPY.double_t ang_y,
    NPY.double_t ang_z,
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center = _default_rotation_center_,# = NPY.array([0.,0.,0.],dtype=NPY.double),
    bint reverse = False):
  """
  apply_euler_angles_and_rotation_center_to_coords(coords,ang_x,ang_y,ang_z,order=[2,1,0],intrinsic=True,rotation_center=[0.,0.,0.],reverse=False)

  Applies the rigid transform corresponding to the Euler angles around the
  rotation center to the provided coords array

  Parameters:
    coords          (np.ndarray[np.double_t]) : Vector of coordinates (shape = (*n*,3))
    ang_x           (double)                  : Rotation angle around the x-axis
    ang_y           (double)                  : Rotation angle around the y-axis
    ang_z           (double)                  : Rotation angle around the z-axis
    order           (np.ndarray[np.int32_t])  : Order of rotations to apply
    intrinsic       (bool)                    : `Axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True, applies the reverse transformation

  Returns:
    Transformed coordinates (`np.ndarray[np.double_t]`, shape = (*n*,3))
  """
  _check_euler_angles_order(order)
  _check_axis_shape(rotation_center,"rotation_center")
  if coords.shape[1] != 3:
    raise AssertionError(f"'coords' argument of invalid shape {NPY.shape(coords)}, expects (n,3)")
  cdef int* order_data = np_to_int_pointer(order)
  cdef int n_samp = coords.shape[0]
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_out = NPY.empty((n_samp,3),dtype=NPY.double)
  PDM_rotation_apply_euler_angles_and_rotation_center(ang_x,
                                                      ang_y,
                                                      ang_z,
                                                      order_data,
                                                      <PDM_bool_t> intrinsic,
                                                      <double*> rotation_center.data,
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
    NPY.ndarray[NPY.int32_t, mode='c', ndim=1] order = _default_order_,# = NPY.array([2,1,0],dtype=NPY.int32),
    bint intrinsic = True,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center = _default_rotation_center_,# = NPY.array([0.,0.,0.],dtype=NPY.double),
    bint reverse = False):
  """
  apply_euler_angles_and_rotation_center_to_vector_field(vector_field,ang_x,ang_y,ang_z,order=[2,1,0],intrinsic=True,rotation_center=[0.,0.,0.],reverse=False)

  Applies the rigid transform corresponding to the Euler angles around the
  rotation center to the provided vector field array

  Parameters:
    vector_field    (np.ndarray[np.double_t]) : Vector (shape = (*n*,3))
    ang_x           (double)                  : Rotation angle around the x-axis
    ang_y           (double)                  : Rotation angle around the y-axis
    ang_z           (double)                  : Rotation angle around the z-axis
    order           (np.ndarray[np.int32_t])  : Order of rotations to apply
    intrinsic       (bool)                    : `Axis conventions <https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations>`_
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True, applies the reverse transformation

  Returns:
    Transformed vector field (`np.ndarray[np.double_t]`, shape = (*n*,3))
  """
  _check_euler_angles_order(order)
  _check_axis_shape(rotation_center,"rotation_center")
  if vector_field.shape[1] != 3:
    raise AssertionError(f"'vector_field' argument of invalid shape {NPY.shape(vector_field)}, expects (n,3)")
  cdef int* order_data = np_to_int_pointer(order)
  cdef double rotation_center_data[3]
  rotation_center_data[:] = [0.,0.,0.]
  cdef int n_samp = vector_field.shape[0]
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_out = NPY.empty((n_samp,3),dtype=NPY.double)
  PDM_rotation_apply_euler_angles_and_rotation_center(ang_x,
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

def apply_axis_angle_and_rotation_center_to_coords(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] coords,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis,
    NPY.double_t angle,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center = _default_rotation_center_,# = NPY.array([0.,0.,0.],dtype=NPY.double),
    bint reverse = False):
  """
  apply_axis_angle_and_rotation_center_to_coords(coords,axis,angle,rotation_center=[0.,0.,0.],reverse=False)

  Applies the rigid transform corresponding to the rotation of the provided angle around the
  provided axis and rotation center to the provided coordinate array

  Parameters:
    coords          (np.ndarray[np.double_t]) : Vector of coordinates (shape = (*n*,3))
    axis            (np.ndarray[np.double_t]) : Rotation axis (shape = (3,))
    angle           (double)                  : Rotation angle (in *radians*)
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True, applies the reverse transformation

  Returns:
    Transformed coordinates (`np.ndarray[np.double_t]`, shape = (*n*,3))
  """
  _check_axis_shape(axis)
  _check_axis_shape(rotation_center,"rotation_center")
  if coords.shape[1] != 3:
    raise AssertionError(f"'coords' argument of invalid shape {NPY.shape(coords)}, expects (n,3)")
  cdef int n_samp = coords.shape[0]
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_out = NPY.empty((n_samp,3),dtype=NPY.double)
  PDM_rotation_apply_axis_angle_and_rotation_center(<double*> axis.data,
                                                    angle,
                                                    <double*> rotation_center.data,
                                                    <PDM_bool_t> reverse,
                                                    <double*> coords.data,
                                                    n_samp,
                                                    <double*> vector_out.data)

  return vector_out

def apply_axis_angle_and_rotation_center_to_vector_field(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_field,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] axis,
    NPY.double_t angle,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center = _default_rotation_center_,# = NPY.array([0.,0.,0.],dtype=NPY.double),
    bint reverse = False):
  """
  apply_axis_angle_and_rotation_center_to_vector_field(vector_field,axis,angle,rotation_center=[0.,0.,0.],reverse=False)

  Applies the rigid transform corresponding to rotation of the provided angle around the
  provided axis and rotation center to the provided vector field array

  Parameters:
    vector_field    (np.ndarray[np.double_t]) : Vector (shape = (*n*,3))
    axis            (np.ndarray[np.double_t]) : Rotation axis (shape = (3,))
    angles          (double)                  : Rotation angle (in *radians*)
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True, applies the reverse transformation

  Returns:
    Transformed vector field (`np.ndarray[np.double_t]`, shape = (*n*,3))
  """
  _check_axis_shape(axis)
  _check_axis_shape(rotation_center)
  if vector_field.shape[1] != 3:
    raise AssertionError(f"'vector_field' argument of invalid shape {NPY.shape(vector_field)}, expects (n,3)")
  cdef int n_samp = vector_field.shape[0]
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_out = NPY.empty((n_samp,3),dtype=NPY.double)
  cdef double rotation_center_data[3]
  rotation_center_data[:] = [0.,0.,0.]
  PDM_rotation_apply_axis_angle_and_rotation_center(<double*> axis.data,
                                                    angle,
                                                    rotation_center_data,
                                                    <PDM_bool_t> reverse,
                                                    <double*> vector_field.data,
                                                    n_samp,
                                                    <double*> vector_out.data)

  return vector_out

def apply_rotation_matrix_and_rotation_center_to_coords(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] coords,
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] rotation_matrix,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center = _default_rotation_center_,# = NPY.array([0.,0.,0.],dtype=NPY.double),
    bint reverse = False):
  """
  apply_rotation_matrix_and_rotation_center_to_coords(coords,rotation_matrix,rotation_center=[0.,0.,0.],reverse=False)

  Applies the rigid transform corresponding to rotation of the provided rotation
  matrix around the provided axis and rotation center to the provided coordinate array

  Parameters:
    coords          (np.ndarray[np.double_t]) : Vector of coordinates (shape = (*n*,3))
    rotation_matrix (np.ndarray[np.double_t]) : 3-by-3 rotation matrix (shape = (3,3))
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True, applies the reverse transformation

  Returns:
    Transformed coordinates (`np.ndarray[np.double_t]`, shape = (*n*,3))
  """
  _check_matrix_shape(rotation_matrix,(3,3))
  _check_axis_shape(rotation_center,"rotation_center")
  if coords.shape[1] != 3:
    raise AssertionError(f"'coords' argument of invalid shape {NPY.shape(coords)}, expects (n,3)")
  cdef int n_samp = coords.shape[0]
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_out = NPY.empty((n_samp,3),dtype=NPY.double)
  PDM_rotation_apply_rotation_matrix_and_rotation_center(<double*> rotation_matrix.data,
                                                         <double*> rotation_center.data,
                                                         <PDM_bool_t> reverse,
                                                         <double*> coords.data,
                                                         n_samp,
                                                         <double*> vector_out.data)
  return vector_out

def apply_homogeneous_matrix_to_coords(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] coords,
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] homogeneous_matrix):
  """
  apply_homogeneous_matrix_to_coords(coords,homogeneous_matrix)

  Applies the rigid transform corresponding to rotation of the provided rotation
  matrix around the provided axis and rotation center to the provided coordinate array

  Parameters:
    coords          (np.ndarray[np.double_t])    : Vector of coordinates (shape = (*n*,3))
    homogeneous_matrix (np.ndarray[np.double_t]) : 4-by-4 rotation matrix (shape = (4,4))

  Returns:
    Transformed coordinates (`np.ndarray[np.double_t]`, shape = (*n*,3))
  """
  if coords.shape[1] != 3:
    raise AssertionError(f"'coords' argument of invalid shape {NPY.shape(coords)}, expects (n,3)")
  _check_matrix_shape(homogeneous_matrix,(4,4),name="homogeneous matrix")
  cdef int n_samp = coords.shape[0]
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_out = NPY.empty((n_samp,3),dtype=NPY.double)
  PDM_rotation_apply_homogeneous_matrix(<double*> homogeneous_matrix.data,
                                        <double*> coords.data,
                                        n_samp,
                                        <double*> vector_out.data)
  return vector_out

def apply_rotation_matrix_and_rotation_center_to_vector_field(
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_field,
    NPY.ndarray[NPY.double_t, mode='c', ndim=2] rotation_matrix,
    NPY.ndarray[NPY.double_t, mode='c', ndim=1] rotation_center = _default_rotation_center_,# = NPY.array([0.,0.,0.],dtype=NPY.double),
    bint reverse = False):
  """
  apply_rotation_matrix_and_rotation_center_to_vector_field(vector_field,rotation_matrix,rotation_center=[0.,0.,0.],reverse=False)

  Applies the rigid transform corresponding to rotation of the provided rotation
  matrix around the provided axis and rotation center to the provided vector field array

  Parameters:
    vector_field    (np.ndarray[np.double_t]) : Vector of coordinates (shape = (*n*,3))
    rotation_matrix (np.ndarray[np.double_t]) : 3-by-3 rotation matrix (shape = (3,3))
    rotation_center (np.ndarray[np.double_t]) : 3D rotation center (shape = (3,))
    reverse         (bool)                    : If True, applies the reverse transformation

  Returns:
    Transformed vector field (`np.ndarray[np.double_t]`, shape = (*n*,3))
  """
  _check_matrix_shape(rotation_matrix,(3,3))
  _check_axis_shape(rotation_center,"rotation_center")
  if vector_field.shape[1] != 3:
    raise AssertionError(f"'vector_field' argument of invalid shape {NPY.shape(vector_field)}, expects (n,3)")
  cdef int n_samp = vector_field.shape[0]
  cdef NPY.ndarray[NPY.double_t, mode='c', ndim=2] vector_out = NPY.empty((n_samp,3),dtype=NPY.double)
  cdef double rotation_center_data[3]
  rotation_center_data[:] = [0.,0.,0.]
  PDM_rotation_apply_rotation_matrix_and_rotation_center(<double*> rotation_matrix.data,
                                                         <double*> rotation_center_data,
                                                         <PDM_bool_t> reverse,
                                                         <double*> vector_field.data,
                                                         n_samp,
                                                         <double*> vector_out.data)
  return vector_out

# endregion
