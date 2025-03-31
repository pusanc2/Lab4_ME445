#!/usr/bin/env python3
import numpy as np
from scipy.linalg import expm
from lab4_header import *

"""
Use 'expm' for matrix exponential.
Angles are in radian, distance are in meters.
"""
def Get_MS():
	# =================== Your code starts here ====================#
	# Fill in the correct values for a1~6 and q1~6, as well as the M matrix
	M = np.array([[0, -1, 0, 390], [0, 0, -1, 401], [1, 0, 0, 215.5], [0, 0, 0, 1]])
	
	w1 = np.array([0, 0, 1])
	w2 = np.array([0, 1, 0])
	w3 = np.array([0, 1, 0])
	w4 = np.array([0, 1, 0])
	w5 = np.array([1, 0, 0])
	w6 = np.array([0, 1, 0])
	q1 = np.array([-150,150,10])
	q2 = np.array([-150,270,162])
	q3 = np.array([94,270,162])
	q4 = np.array([307,177,162])
	q5 = np.array([307,260,162])	
	q6 = np.array([390,260,162])
	v1 = np.cross(-w1, q1)
	v2 = np.cross(-w2, q2)
	v3 = np.cross(-w3, q3)
	v4 = np.cross(-w4, q4)
	v5 = np.cross(-w5, q5)
	v6 = np.cross(-w6, q6)
	S1 = np.array([[0, -w1[2], w1[1], v1[0]], [w1[2], 0, -w1[0], v1[1]], [-w1[1], w1[0], 0, v1[2]], [0, 0, 0, 0]])
	S2 = np.array([[0, -w2[2], w2[1], v2[0]], [w2[2], 0, -w2[0], v2[1]], [-w2[1], w2[0], 0, v2[2]], [0, 0, 0, 0]])
	S3 = np.array([[0, -w3[2], w3[1], v3[0]], [w3[2], 0, -w3[0], v3[1]], [-w3[1], w3[0], 0, v3[2]], [0, 0, 0, 0]])
	S4 = np.array([[0, -w4[2], w4[1], v4[0]], [w4[2], 0, -w4[0], v4[1]], [-w4[1], w4[0], 0, v4[2]], [0, 0, 0, 0]])
	S5 = np.array([[0, -w5[2], w5[1], v5[0]], [w5[2], 0, -w5[0], v5[1]], [-w5[1], w5[0], 0, v5[2]], [0, 0, 0, 0]])
	S6 = np.array([[0, -w6[2], w6[1], v6[0]], [w6[2], 0, -w6[0], v6[1]], [-w6[1], w6[0], 0, v6[2]], [0, 0, 0, 0]])
	
	L1 = 152
	L2 = 120
	L3 = 244
	L4 = 93
	L5 = 213
	L6 = 83
	L7 = 83
	L8 = 82
	L9 = 53.5
	print(M, "\n")

	# ==============================================================#
	return M, S1, S2, S3, S4, S5, S6


"""
Function that calculates encoder numbers for each motor
"""
def lab_fk(theta1, theta2, theta3, theta4, theta5, theta6):

	# Initialize the return_value
	return_value = [None, None, None, None, None, None]

	print("Foward kinematics calculated:\n")

	# =================== Your code starts here ====================#
	theta = np.array([theta1,theta2,theta3,theta4,theta5,theta6])
	T = np.eye(4)

	M, S1, S2, S3, S4, S5, S6= Get_MS()
	T = expm(S1*theta1)@expm(S2*theta2)@expm(S3*theta3)@expm(S4*theta4)@expm(S5*theta5)@expm(S6*theta6)@M



	# ==============================================================#
	print(str(T) + "\n")

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
	L1 = 152
	L2 = 120
	L3 = 244
	L4 = 93
	L5 = 213
	L6 = 83
	L7 = 83
	L8 = 82
	L9 = 53.5
	print(xWgrip)
	xWgrip = xWgrip*1000 + 150
	yWgrip = yWgrip*1000 - 150
	zWgrip = zWgrip*1000 - 10
	yaw_WgripDegree = np.deg2rad(yaw_WgripDegree)
	xcen =  xWgrip - (L9*np.cos(yaw_WgripDegree))
	ycen = yWgrip - (L9*np.sin(yaw_WgripDegree))
	zcen  = zWgrip
	
	a = np.sqrt(xcen**2 + ycen**2)
	b = L2-L4+L6
	print(a)
	print(b)
	stheta = np.arcsin(b/a)
	total_theta = np.arctan(ycen/xcen)
	theta1 = total_theta - stheta
	# print(theta1)
	theta6 = (np.pi/2)-yaw_WgripDegree+theta1
 #####################################################################
	x3end = xcen - np.cos(theta1)*L7 + np.sin(theta1)*110 # 83+27
	y3end = ycen - np.sin(theta1)*L7 + np.cos(theta1)*110
	z3end = zcen + 59 + L8
	
 
	r = np.sqrt(x3end**2 + y3end**2)
	 
	delta_z = z3end - L1
	d = np.sqrt(r**2 + delta_z**2)
	v = np.arccos((L5**2 + L3**2 - d**2) / (2 * L5 * L3))
	theta3 = np.pi - v
	theta2_1 = np.arccos(r / d)
	theta2_2 = np.arccos((L3**2 + d**2 - L5**2) / (2 * L3 * d))
	
	theta2 = -(theta2_1 + theta2_2)
 
	
	theta4 = -(theta2 + theta3)
	theta5 = -np.pi/2
	print("theta1: ", theta1 * 180 / np.pi, "\ntheta2: ", theta2 * 180 / np.pi, "\ntheta3: ", theta3 * 180 / np.pi, "\ntheta4: ", theta4 * 180 / np.pi, "\ntheta5: ", theta5 * 180 / np.pi, "\ntheta6: ", theta6 * 180 / np.pi)
	# ==============================================================#
	return lab_fk(theta1, theta2, theta3, theta4, theta5, theta6)
    
    
    
    
    
    