#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np

import Pypdm.Pypdm as PDM

comm = MPI.COMM_WORLD

i_rank = comm.rank
n_rank = comm.size

vector = np.array([
    [1.,0.,0.],
    [0.,1.,0.],
    [0.,0.,1.],
    [1.,1.,1.],
],dtype=np.float64)
deg2rad = np.pi/180.
rotation_center = np.array([1.,2.,3.],dtype=np.float64)

ang_x = 5. *deg2rad
ang_y = 10.*deg2rad
ang_z = 15.*deg2rad

#region conversions ------------------------------------------------------------

print(PDM.axis_angle_to_euler_angles.__doc__)
print(PDM.axis_angle_to_rotation_matrix.__doc__)
print(PDM.euler_angles_to_axis_angle.__doc__)
print(PDM.euler_angles_to_euler_angles.__doc__)
print(PDM.euler_angles_to_rotation_matrix.__doc__)
print(PDM.rotation_matrix_to_axis_angle.__doc__)
print(PDM.rotation_matrix_to_euler_angles.__doc__)
print(PDM.two_vectors_to_axis_angle.__doc__)
print(PDM.two_vectors_to_euler_angles.__doc__)
print(PDM.two_vectors_to_rotation_matrix.__doc__)
exit()
# not implemented yet
# PDM.two_vectors_to_axis_angle
# PDM.two_vectors_to_euler_angles
# PDM.two_vectors_to_rotation_matrix

#region euler angles -----------------------------------------------------------


## apply_euler_angles_and_rotation_center_to_coords ---

vector_out = PDM.apply_euler_angles_and_rotation_center_to_coords(vector,
                                                        ang_x,
                                                        ang_y,
                                                        ang_z,
                                                        rotation_center=rotation_center)
print(vector_out)
#[[ 0.91747916  0.18590612 -0.11484409]
# [-0.2769875   0.89718638  0.14463574]
# [ 0.15587848 -0.10839444  1.03986435]
# [ 0.8639143   1.11265983  0.95204782]]
vector_out_out = PDM.apply_euler_angles_and_rotation_center_to_coords(vector_out,
                                                        ang_x,
                                                        ang_y,
                                                        ang_z,
                                                        rotation_center=rotation_center,
                                                        reverse=True)
print(vector_out_out)
# [[ 1.00000000e+00 -2.42861287e-16 -2.42861287e-16]
#  [-4.85722573e-17  1.00000000e+00 -2.28983499e-16]
#  [-4.85722573e-17 -2.35922393e-16  1.00000000e+00]
#  [ 1.00000000e+00  1.00000000e+00  1.00000000e+00]]
print(np.allclose(vector_out_out,vector))
# True

## apply_euler_angles_and_rotation_center_to_vector_field ---
vector_out = PDM.apply_euler_angles_and_rotation_center_to_vector_field(vector,
                                                        ang_x,
                                                        ang_y,
                                                        ang_z,
                                                        rotation_center=rotation_center)
print(vector_out)
#[[ 0.95125124  0.254887   -0.17364818]
# [-0.24321542  0.96616727  0.08583165]
# [ 0.18965056 -0.03941355  0.98106026]
# [ 0.89768638  1.18164072  0.89324374]]
vector_out_out = PDM.apply_euler_angles_and_rotation_center_to_vector_field(vector_out,
                                                        ang_x,
                                                        ang_y,
                                                        ang_z,
                                                        rotation_center=rotation_center,
                                                        reverse=True)
print(vector_out_out)
# [[ 1.00000000e+00 -3.61989592e-17 -1.39611063e-17]
#  [-3.61989592e-17  1.00000000e+00  4.64042853e-18]
#  [-1.39611063e-17  4.64042853e-18  1.00000000e+00]
#  [ 1.00000000e+00  1.00000000e+00  1.00000000e+00]]
print(np.allclose(vector_out_out,vector))
# True

#endregion
#region axis angle -------------------------------------------------------------

axis = np.array([1.,2.,3.])
angle = 12.3456*deg2rad

vector_out = PDM.apply_axis_angle_and_rotation_center_to_coords(vector,
                                                            axis,
                                                            angle,
                                                            rotation_center=rotation_center)
print(vector_out)
# [[ 0.97852745  0.17473118 -0.10932994]
#  [-0.16812424  0.98348266  0.06705298]
#  [ 0.11924034 -0.04723216  0.99174133]
#  [ 0.92964356  1.11098167  0.94946437]]
vector_out_out = PDM.apply_axis_angle_and_rotation_center_to_coords(vector_out,
                                                            axis,
                                                            angle,
                                                            rotation_center=rotation_center,
                                                            reverse=True)
print(vector_out_out)
# [[1.00000000e+00 2.30559585e-16 2.56660610e-19]
#  [1.19537283e-16 1.00000000e+00 1.58216112e-17]
#  [9.74011753e-17 2.10110640e-16 1.00000000e+00]
#  [1.00000000e+00 1.00000000e+00 1.00000000e+00]]

print(np.allclose(vector_out_out,vector))

#endregion
#region rotation matrix --------------------------------------------------------

rotation_matrix = PDM.euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z)

vector_out = PDM.apply_rotation_matrix_and_rotation_center_to_coords(vector,
                                                            rotation_matrix,
                                                            rotation_center=rotation_center)
print(vector_out)
# [[ 0.91747916  0.18590612 -0.11484409]
#  [-0.2769875   0.89718638  0.14463574]
#  [ 0.15587848 -0.10839444  1.03986435]
#  [ 0.8639143   1.11265983  0.95204782]]

vector_out_out = PDM.apply_rotation_matrix_and_rotation_center_to_coords(vector_out,
                                                            rotation_matrix,
                                                            rotation_center=rotation_center,
                                                            reverse=True)
print(vector_out_out)
# [[ 1.00000000e+00 -2.42861287e-16 -2.42861287e-16]
#  [-4.85722573e-17  1.00000000e+00 -2.28983499e-16]
#  [-4.85722573e-17 -2.35922393e-16  1.00000000e+00]
#  [ 1.00000000e+00  1.00000000e+00  1.00000000e+00]]
print(np.allclose(vector_out_out,vector))
#endregion

#region apply rotation to all the components of a pytree -----------------------

# coordinates
# vector field
# other invariant variables
# periodic -> convert to axis-angle -> rotate axis -> convert back to euler angles


print("[{}] -- End".format(i_rank))