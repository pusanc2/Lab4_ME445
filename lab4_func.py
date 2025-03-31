#!/usr/bin/env python3
import numpy as np
from scipy.linalg import expm
from lab4_header import *

# Test 1 (100, 100, 150, 90):
# x = 99 mm
# y = 103 mm
# z = 145 mm
#
# Test 2 (100, 500, 100, 90):
# x = 98 mm
# y = 504 mm
# z = 96 mm

# helper functions for lab_fk
def VecTose3(V):
	return np.r_[np.c_[VecToso3([V[0], V[1], V[2]]), [V[3], V[4], V[5]]], np.zeros((1,4))]

def VecToso3(omg):
	return np.array([[0, -omg[2], omg[1]],[omg[2], 0, -omg[0]], [-omg[1], omg[0], 0]])



"""
Use 'expm' for matrix exponential.
Angles are in radian, distance are in meters.
"""
def Get_MS():
	# =================== Your code starts here ====================#
	# Fill in the correct values for a1~6 and q1~6, as well as the M matrix
 
	M = np.array([[0, -1, 0, -150 + 244 + 213 + 83],
               	  [0, 0, -1, 150 + 120 - 93 + 83 + 82 + 59],
                  [1, 0,  0, 215.5],
                  [0, 0,  0, 1],])
	
 
	q1 = np.array([-150, 				  150, 				   10])
	q2 = np.array([-150, 				  150 + 120, 		   10 + 152])
	q3 = np.array([-150 + 244,			  150 + 120, 		   10 + 152])
	q4 = np.array([-150 + 244 + 213,	  150 + 120 - 93, 	   10 + 152])
	q5 = np.array([-150 + 244 + 213,      150 + 120 - 93 + 83, 10 + 152])
	q6 = np.array([-150 + 244 + 213 + 83, 150 + 120 - 93 + 83, 10 + 152])

	w1 = np.array([0,0,1])
	w2 = np.array([0,1,0])
	w3 = np.array([0,1,0])
	w4 = np.array([0,1,0])
	w5 = np.array([1,0,0])
	w6 = np.array([0,1,0])
	
	S1 = np.block([w1, np.cross(-w1,q1)])
	S2 = np.block([w2, np.cross(-w2,q2)])
	S3 = np.block([w3, np.cross(-w3,q3)])
	S4 = np.block([w4, np.cross(-w4,q4)])
	S5 = np.block([w5, np.cross(-w5,q5)])
	S6 = np.block([w6, np.cross(-w6,q6)])
 
	S = np.block([[S1], [S2], [S3], [S4], [S5], [S6]])
 


	# ==============================================================#
	return M, S


"""
Function that calculates encoder numbers for each motor
"""
def lab_fk(theta1, theta2, theta3, theta4, theta5, theta6):

	# Initialize the return_value
	return_value = [None, None, None, None, None, None]

	print("Thetas calculated:\n")
 

	# =================== Your code starts here ====================#
	theta = np.degrees(np.array([theta1,theta2,theta3,theta4,theta5,theta6]))
 
	print(str(theta) + '\n')
 
	T = np.eye(4)
	M, S = Get_MS()
 
	S1_mat = VecTose3(S[0])
	S2_mat = VecTose3(S[1])
	S3_mat = VecTose3(S[2])
	S4_mat = VecTose3(S[3])
	S5_mat = VecTose3(S[4])
	S6_mat = VecTose3(S[5])
 
	T = expm(S1_mat * theta1) @ expm(S2_mat * theta2) @ expm(S3_mat * theta3) @ expm(S4_mat * theta4) @ expm(S5_mat * theta5) @ expm(S6_mat * theta6) @ M
	
	print("Foward kinematics calculated:\n")
	print(str(T) + "\n")


	# ==============================================================#

	return_value[0] = theta1 + PI
	return_value[1] = theta2
	return_value[2] = theta3
	return_value[3] = theta4 - (0.5*PI)
	return_value[4] = theta5
	return_value[5] = theta6

	return return_value


"""
Function that calculates an elbow up Inverse Kinematic solution for the UR3
"""
def lab_invk(xWgrip, yWgrip, zWgrip, yaw_WgripDegree):
	# =================== Your code starts here ====================#
	
	# convert World coordinates to Base coordinates (frame 0)
	x_grip = xWgrip + 150
	y_grip = yWgrip - 150
	z_grip = zWgrip - 10
	yaw_WgripRadians = np.radians(yaw_WgripDegree)
	
	# arm lengths
	L1 = 152
	L2 = 120
	L3 = 244
	L4 = 93
	L5 = 213
	L6 = 83
	L7 = 83
	L8 = 82
	L9 = 53.5
	L10 = 59
 
	# find wrist's center point
	x_cen = x_grip - L9*np.cos(yaw_WgripRadians)
	y_cen =	y_grip - L9*np.sin(yaw_WgripRadians)
	z_cen = z_grip
	
	# find theta1 (waist angle)
	if (((L2-L4+L6)/np.sqrt(x_cen**2 + y_cen**2)) > 1):
		print("WARNING: INVALID POSITION")
 
	theta_small = np.arcsin((L2-L4+L6)/np.sqrt(x_cen**2 + y_cen**2))
	theta_big = np.arctan2(y_cen, x_cen)
	theta1 = theta_big - theta_small
 
	
	# print("x_cen, y_cen, z_cen:", x_cen, y_cen, z_cen)
	# print("theta_small:", np.degrees(theta_small), "theta_big:", np.degrees(theta_big))
	
 
	# find theta6
	theta6 = np.pi/2 + theta1 - yaw_WgripRadians
 
	# find projected endpoint
	x_3end = x_cen - L7*np.cos(theta1) + (L6+27)*np.sin(theta1)
	y_3end = y_cen - L7*np.sin(theta1) - (L6+27)*np.cos(theta1)
	z_3end = z_cen + L10 + L8 
 
	# find theta2, theta3, theta4
	hypotenuse_xy = np.sqrt((x_3end)**2 + (y_3end)**2)	
	
	d = np.sqrt((z_3end - L1)**2 + (hypotenuse_xy)**2)
	thetau = np.arccos((L3**2 + L5**2 - d**2)/(2*L3*L5))
	
	theta3 = np.pi - thetau
 
	theta2_big = np.arccos((L3**2 + d**2 - L5**2)/(2*L3*d))
	theta2_small = np.arctan2((z_3end - L1), (hypotenuse_xy))
	theta2 = -(theta2_big + theta2_small)
 
	theta4 = -(theta3 + theta2)
 
 
	theta5 = -np.pi/2
	# ==============================================================#
	return lab_fk(theta1, 
               	theta2, 
                theta3, 
                theta4, 
                theta5, 
                theta6)
	# return lab_fk(np.degrees(theta1), 
    #            	np.degrees(theta2), 
    #             np.degrees(theta3), 
    #             np.degrees(theta4), 
    #             np.degrees(theta5), 
    #             np.degrees(theta6))
