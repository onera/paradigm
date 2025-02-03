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

class Test_axis_angle_to_euler_angles:
    @staticmethod
    def test_arguments():
        axes = [
            (np.array([1.,0.,0.]),True),
            (np.array([2.,0.,0.]),True),
            (np.array([1.,1.,0.]),True),
            (np.array([0.,0.,1.]),True),
            (np.array([0.,0.,1.,0.]),False),
            (np.array([[0.,0.],[1.,0.]]),False),
            (np.array([[0.,0.],[1.,0.]]),False),
            ("toto",False),
        ]
        angles = [
            (1.,True),
            ("toto",False)
        ]
        orders = [
            (np.array([2,1,0],dtype=np.int32),True),
            (np.array([1,2,0],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int64),False),
            ("toto",False),
            (np.array(0.14),False),
            (np.array([0.,1.,0.,0.]),False),
            (np.array([0,1,2,3]),False),
            (np.array([1,1,0]),False),
        ]
        intrinsic = [
            (True,True),
            (False,True),
            ("toto",False)
        ]
        for axis,axis_valid in axes:
            for angle,angle_valid in angles:
                if axis_valid and angle_valid:
                    out = PDM.axis_angle_to_euler_angles(axis,angle)  
                else:
                    try:
                        out = PDM.axis_angle_to_euler_angles(axis,angle)  
                    except:
                        pass

        for order,order_valid in orders:
            for intr,intrinsic_valid in intrinsic:
                if order_valid and intrinsic_valid:
                    out = PDM.axis_angle_to_euler_angles(np.array([1.,1.,2.]),12,order,intr)
                else:
                    try:
                        out = PDM.axis_angle_to_euler_angles(np.array([1.,1.,2.]),12,order,intr)
                    except:
                        pass

    @staticmethod
    def test_compute():
        axis = np.array([2.,0.,0.])
        angle = 45*deg2rad
        out = PDM.axis_angle_to_euler_angles(axis,angle)  
        expec = (angle,0.,0.)
        import itertools
        for intr in [True,False]:
            for order in itertools.permutations(range(3),3):
                order = np.array(order,dtype=np.int32)
                out = PDM.axis_angle_to_euler_angles(axis,angle,order=order,intrinsic=intr)  
                for a,ex in zip(out,expec):
                    assert abs(a-ex)<1e-12

class Test_axis_angle_to_rotation_matrix:

    @staticmethod
    def test_arguments():
        axes = [
            (np.array([1.,0.,0.]),True),
            (np.array([2.,0.,0.]),True),
            (np.array([1.,1.,0.]),True),
            (np.array([0.,0.,1.]),True),
            (np.array([0.,0.,1.,0.]),False),
            (np.array([[0.,0.],[1.,0.]]),False),
            (np.array([[0.,0.],[1.,0.]]),False),
            ("toto",False),
        ]
        angles = [
            (1.,True),
            ("toto",False)
        ]
        for axis,axis_valid in axes:
            for angle,angle_valid in angles:
                if axis_valid and angle_valid:
                    out = PDM.axis_angle_to_rotation_matrix(axis,angle)  
                else:
                    try:
                        out = PDM.axis_angle_to_rotation_matrix(axis,angle)  
                    except:
                        pass

    @staticmethod
    def test_compute():
        axis = np.array([2.,1.,5.])
        angle = 12*deg2rad
        out = PDM.axis_angle_to_rotation_matrix(axis,angle)
        expec = np.array([
            [ 0.98106125, -0.18833971,  0.04524344],
            [ 0.19125337,  0.97887601, -0.07227655],
            [-0.03067517,  0.07956068,  0.99635793]])
        assert np.allclose(out,expec)

class Test_euler_angles_to_axis_angle:
    @staticmethod
    def test_arguments():
        angles = [
            ((1.,0.,0.),True),
            (("toto",0.,1.),False),
            ((1.,"toto",0.),False),
            ((1.,0.,"toto"),False),
        ]
        orders = [
            (np.array([2,1,0],dtype=np.int32),True),
            (np.array([1,2,0],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int64),False),
            ("toto",False),
            (np.array(0.14),False),
            (np.array([0.,1.,0.,0.]),False),
            (np.array([0,1,2,3]),False),
            (np.array([1,1,0]),False),
        ]
        intrinsic = [
            (True,True),
            (False,True),
            ("toto",False)
        ]
        for angl,angle_valid in angles:
            for order,order_valid in orders:
                for intr,intrinsic_valid in intrinsic:
                    if angle_valid and order_valid and intrinsic_valid:
                        out = PDM.euler_angles_to_axis_angle(*angl,order,intr)  
                    else:
                        try:
                            out = PDM.euler_angles_to_axis_angle(*angl,order,intr)
                        except:
                            pass
    @staticmethod
    def test_compute():
        ang_x = 12.*deg2rad
        ang_y = 20.*deg2rad
        ang_z = 4. *deg2rad
        axis,angle = PDM.euler_angles_to_axis_angle(ang_x,ang_y,ang_z)
        
        expec_axis = np.array([0.4801993,  0.8735463,  0.07953301])
        expec_angle = 0.4061631775316828
        assert np.allclose(axis,expec_axis)
        assert np.allclose(angle,expec_angle)

        rotation_matrix = PDM.euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z)
        axis2,angle2 = PDM.rotation_matrix_to_axis_angle(rotation_matrix)
        assert np.allclose(axis2,expec_axis)
        assert np.allclose(angle2,expec_angle)


        axis,angle = PDM.euler_angles_to_axis_angle(ang_x,ang_y,ang_z,intrinsic=False)  
        expec_axis = np.array([0.52422041, 0.81348871, 0.2518513])
        expec_angle = 0.41854053477833275
        assert np.allclose(axis,expec_axis)
        assert np.allclose(angle,expec_angle)
        
        rotation_matrix = PDM.euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z,intrinsic=False)
        axis2,angle2 = PDM.rotation_matrix_to_axis_angle(rotation_matrix)
        assert np.allclose(axis2,expec_axis)
        assert np.allclose(angle2,expec_angle)
        
class Test_euler_angles_to_euler_angles:
    @staticmethod
    def test_arguments():
        angles = [
            ((1.,0.,0.),True),
            (("toto",0.,1.),False),
            ((1.,"toto",0.),False),
            ((1.,0.,"toto"),False),
        ]
        orders = [
            (np.array([2,1,0],dtype=np.int32),True),
            (np.array([1,2,0],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int64),False),
            ("toto",False),
            (np.array(0.14),False),
            (np.array([0.,1.,0.,0.]),False),
            (np.array([0,1,2,3]),False),
            (np.array([1,1,0]),False),
        ]
        intrinsic = [
            (True,True),
            (False,True),
            ("toto",False)
        ]
        for angl,angle_valid in angles:
            for order,order_valid in orders:
                for intr,intrinsic_valid in intrinsic:
                    for o_order,o_order_valid in orders:
                        for o_intr,o_intrinsic_valid in intrinsic:
                            if angle_valid and order_valid and intrinsic_valid and o_order_valid and o_intrinsic_valid:
                                out = PDM.euler_angles_to_euler_angles(*angl,order,intr,o_order,o_intr)
                            else:
                                try:
                                    out = PDM.euler_angles_to_euler_angles(*angl,order,intr,o_order,o_intr)
                                except:
                                    pass

    @staticmethod
    def test_compute():
        ang_x = 12.*deg2rad
        ang_y = 20.*deg2rad
        ang_z = 4. *deg2rad
        order = np.array([2,1,0],dtype=np.int32)
        o_order = np.array([0,1,2],dtype=np.int32)
        out = PDM.euler_angles_to_euler_angles(ang_x,ang_y,ang_z,order,True,o_order,True)
        expec = (0.19764331809693458, 0.3556869840140926, -0.002885236490910617)
        for a,e in zip(out,expec):
            assert abs(a-e)<1e-12
        out = PDM.euler_angles_to_euler_angles(ang_x,ang_y,ang_z,order,True,o_order,False)
        assert abs(out[0]-ang_x)<1e-12
        assert abs(out[1]-ang_y)<1e-12
        assert abs(out[2]-ang_z)<1e-12
        
class Test_euler_angles_to_rotation_matrix:
    @staticmethod
    def test_arguments():
        angles = [
            ((1.,0.,0.),True),
            (("toto",0.,1.),False),
            ((1.,"toto",0.),False),
            ((1.,0.,"toto"),False),
        ]
        orders = [
            (np.array([2,1,0],dtype=np.int32),True),
            (np.array([1,2,0],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int64),False),
            ("toto",False),
            (np.array(0.14),False),
            (np.array([0.,1.,0.,0.]),False),
            (np.array([0,1,2,3]),False),
            (np.array([1,1,0]),False),
        ]
        intrinsic = [
            (True,True),
            (False,True),
            ("toto",False)
        ]
        for angl,angle_valid in angles:
            for order,order_valid in orders:
                for intr,intrinsic_valid in intrinsic:
                    if angle_valid and order_valid and intrinsic_valid:
                        out = PDM.euler_angles_to_rotation_matrix(*angl,order,intr)  
                    else:
                        try:
                            out = PDM.euler_angles_to_rotation_matrix(*angl,order,intr)
                        except:
                            pass
    @staticmethod
    def test_compute():
        orders = [
            np.array([2,1,0],dtype=np.int32),
            np.array([1,2,0],dtype=np.int32),
            np.array([2,0,1],dtype=np.int32),
        ]
        for d in range(3):
            ang_x = (12. if d == 0 else 0.)*deg2rad
            ang_y = (12. if d == 1 else 0.)*deg2rad
            ang_z = (12. if d == 2 else 0.)*deg2rad
            for order in orders:
                for intr in [True,False]:
                    out = PDM.euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z,order,intr)
                    for i in range(3):
                        assert out[i,i] == 1. if d==i else np.cos(12.)

        ang_x = 12.*deg2rad
        ang_y = 4. *deg2rad
        ang_z = 20.*deg2rad
        out = PDM.euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z)
        print(out)
        expec = np.array([
            [ 0.93740358, -0.32091765,  0.13522721],
            [ 0.341187  ,  0.92411846, -0.17203632],
            [-0.06975647,  0.20740523,  0.97576488],
        ])
        assert np.allclose(out,expec)

class Test_rotation_matrix_to_axis_angle:
    @staticmethod
    def test_arguments():
        rot_mat = [
            (np.eye(3,3,dtype=np.float64),True),
            (np.eye(4,4,dtype=np.float64),False),
            (np.eye(3,3,dtype=np.int64),False),
            (np.array([1.]),False),
            (np.array([[1],[2]]),False),
        ]
        for mat,mat_valid in rot_mat:
            if mat_valid:
                out = PDM.rotation_matrix_to_axis_angle(mat)
            else:
                try:
                    out = PDM.rotation_matrix_to_axis_angle(mat)
                except:
                    pass

    @staticmethod
    def test_compute():
        axis = np.array([2.,1.,5.])
        angle = 12*deg2rad
        mat = PDM.axis_angle_to_rotation_matrix(axis,angle)
        axis2,angle2 = PDM.rotation_matrix_to_axis_angle(mat)
        axis /=np.linalg.norm(axis)
        assert np.allclose(axis,axis2)
        assert np.allclose(angle,angle2)
        
class Test_rotation_matrix_to_euler_angles:
    @staticmethod
    def test_arguments():
        rot_mat = [
            (np.eye(3,3,dtype=np.float64),True),
            (np.eye(4,4,dtype=np.float64),False),
            (np.eye(3,3,dtype=np.int64),False),
            (np.array([1.]),False),
            (np.array([[1],[2]]),False),
        ]
        orders = [
            (np.array([2,1,0],dtype=np.int32),True),
            (np.array([1,2,0],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int64),False),
            ("toto",False),
            (np.array(0.14),False),
            (np.array([0.,1.,0.,0.]),False),
            (np.array([0,1,2,3]),False),
            (np.array([1,1,0]),False),
        ]
        intrinsic = [
            (True,True),
            (False,True),
            ("toto",False)
        ]
        for mat,mat_valid in rot_mat:
            for order,order_valid in orders:
                for intr,intrinsic_valid in intrinsic:
                    if mat_valid and order_valid and intrinsic_valid:
                        out = PDM.rotation_matrix_to_euler_angles(mat,order,intr)
                    else:
                        try:
                            out = PDM.rotation_matrix_to_euler_angles(mat,order,intr)
                        except:
                            pass

    @staticmethod
    def test_compute():
        orders = [
            np.array([2,1,0],dtype=np.int32),
            np.array([1,2,0],dtype=np.int32),
            np.array([2,0,1],dtype=np.int32),
        ]
        ang_x = 12.*deg2rad
        ang_y = 4. *deg2rad
        ang_z = 20.*deg2rad
        for order in orders:
            for intr in [True,False]:
                mat = PDM.euler_angles_to_rotation_matrix(ang_x,ang_y,ang_z,order,intr)
                ang_x2,ang_y2,ang_z2 = PDM.rotation_matrix_to_euler_angles(mat,order,intr)
                assert abs(ang_x-ang_x2)<1e-12
                assert abs(ang_y-ang_y2)<1e-12
                assert abs(ang_z-ang_z2)<1e-12
        
class Test_two_vectors_to_axis_angle:
    @staticmethod
    def test_arguments():
        vectors = [
            (np.array([1.,0.,0.]),True),
            (np.array([1.,0.,0.,1.]),False),
            (np.array([1.,0.]),False),
            (np.array([[1.,0.]]),False),
            (np.array([[1.,0.,0.],[2.,0.,0.]]),False),
        ]
        for vec1,vec1_valid in vectors:
            for vec2,vec2_valid in vectors:
                if vec1_valid and vec2_valid:
                    out = PDM.two_vectors_to_axis_angle(vec1,vec2)
                else:
                    try:
                        out = PDM.two_vectors_to_axis_angle(vec1,vec2)
                    except:
                        pass

    @staticmethod
    def test_compute():
        vec1 = np.array([1.,0.,0.])
        vec2 = np.array([1.,1.,0.])
        axis,angle = PDM.two_vectors_to_axis_angle(vec1,vec2)
        print(axis,angle)

class Test_two_vectors_to_euler_angles:
    @staticmethod
    def test_arguments():
        vectors = [
            (np.array([1.,0.,0.]),True),
            (np.array([1.,0.,0.,1.]),False),
            (np.array([1.,0.]),False),
            (np.array([[1.,0.]]),False),
            (np.array([[1.,0.,0.],[2.,0.,0.]]),False),
        ]
        orders = [
            (np.array([2,1,0],dtype=np.int32),True),
            (np.array([1,2,0],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int32),True),
            (np.array([2,0,1],dtype=np.int64),False),
            ("toto",False),
            (np.array(0.14),False),
            (np.array([0.,1.,0.,0.]),False),
            (np.array([0,1,2,3]),False),
            (np.array([1,1,0]),False),
        ]
        intrinsic = [
            (True,True),
            (False,True),
            ("toto",False)
        ]
        for vec1,vec1_valid in vectors:
            for vec2,vec2_valid in vectors:
                for order,order_valid in orders:
                    for intr,intrinsic_valid in intrinsic:
                        if vec1_valid and vec2_valid and order_valid and intrinsic_valid:
                            out = PDM.two_vectors_to_euler_angles(vec1,vec2,order,intr)
                        else:
                            try:
                                out = PDM.two_vectors_to_euler_angles(vec1,vec2,order,intr)
                            except:
                                pass
        
    @staticmethod
    def test_compute():
        vec1 = np.array([1.,0.,0.])
        vec2 = np.array([1.,1.,0.])
        out = PDM.two_vectors_to_euler_angles(vec1,vec2)
        print(out)
        
class Test_two_vectors_to_rotation_matrix:
    @staticmethod
    def test_arguments():
        vectors = [
            (np.array([1.,0.,0.]),True),
            (np.array([1.,0.,0.,1.]),False),
            (np.array([1.,0.]),False),
            (np.array([[1.,0.]]),False),
            (np.array([[1.,0.,0.],[2.,0.,0.]]),False),
        ]
        for vec1,vec1_valid in vectors:
            for vec2,vec2_valid in vectors:
                if vec1_valid and vec2_valid:
                    out = PDM.two_vectors_to_rotation_matrix(vec1,vec2)
                else:
                    try:
                        out = PDM.two_vectors_to_rotation_matrix(vec1,vec2)
                    except:
                        pass

    @staticmethod
    def test_compute():
        vec1 = np.array([1.,0.,0.])
        vec2 = np.array([1.,1.,0.])
        out = PDM.two_vectors_to_rotation_matrix(vec1,vec2)
        print(out)

Test_axis_angle_to_euler_angles.test_arguments()
Test_axis_angle_to_euler_angles.test_compute()
Test_axis_angle_to_rotation_matrix.test_arguments()
Test_axis_angle_to_rotation_matrix.test_compute()
Test_euler_angles_to_axis_angle.test_arguments()
Test_euler_angles_to_axis_angle.test_compute()
Test_euler_angles_to_euler_angles.test_arguments()
Test_euler_angles_to_euler_angles.test_compute()
Test_euler_angles_to_rotation_matrix.test_arguments()
Test_euler_angles_to_rotation_matrix.test_compute()
Test_rotation_matrix_to_axis_angle.test_arguments()
Test_rotation_matrix_to_axis_angle.test_compute()
Test_rotation_matrix_to_euler_angles.test_arguments()
Test_rotation_matrix_to_euler_angles.test_compute()
Test_two_vectors_to_axis_angle.test_arguments()
Test_two_vectors_to_axis_angle.test_compute()
Test_two_vectors_to_euler_angles.test_arguments()
Test_two_vectors_to_euler_angles.test_compute()
Test_two_vectors_to_rotation_matrix.test_arguments()
Test_two_vectors_to_rotation_matrix.test_compute()

print(PDM.axis_angle_and_rotation_center_to_homogeneous_matrix)
print(PDM.euler_angles_and_rotation_center_to_homogeneous_matrix)
print(PDM.rotation_matrix_and_rotation_center_to_homogeneous_matrix)
print(PDM.periodic_t_info_to_homogeneous_matrix)
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