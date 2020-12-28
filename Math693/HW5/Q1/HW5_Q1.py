#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Program written for Homework 5, MATH 693B. Problem outline from
# Strikwerda 7.3.3

# Program solves a two-dimensional differential equation using the
# Peaceman-Rachford ADI algorithm. Spatial derivatives calculated using the
# backward-space implicit method.

# Author: Jon Parsons
# Date: 4-13-20

################################################################################
#Functions
################################################################################
# Function containing the initial condition and exact solution of ODE u
def u_exact(x,y):
    # x - point of interest in x
    # y - point of interest in y

    if abs(x) <= 0.5 and abs(y) <= 0.5:
        sln = (1.0 - 2.0*abs(x))*(1.0 - 2.0*abs(y))
    else:
        sln = 0.0

    return sln

################################################################################
# Function that fills the A1 and A2 arrays of the Peaceman-Rachford method
# Takes those and finalizes the 2 matrices for the final computations
def mat_fill(h,k,n):
    # h - dx
    # k - dt
    # n - number of elements

    np.set_printoptions(linewidth=150)

    A1 = np.zeros((round(n*n),round(n*n))) # Matrix containing FE in dx
    A2 = np.zeros((round(n*n),round(n*n))) # Matrix containing FE in dy
    A_sub = np.zeros((round(n),round(n))) # Sub-matrix of A1
    B_sub = np.zeros((round(n),round(n))) # Sub-matrix of A2 (first and last)
    I = np.identity(n*n) # Identity matrix of appropriate size

    # Fill sub-matrices
    for i in range(round(n)):
        B_sub[i,i] = 2*k/(2*h)
        if i != 0 and i != (n-1):
            A_sub[i,i] = k/(2*h)
            A_sub[i,i-1] = -k/(2*h)

    # Fill full matrices
    for i in range(round(n-1)):
        A1[i*n:i*n+n,i*n:i*n+n] = A_sub
        if i > 0 and i < (n-1):
            A2[i*n:i*n+n,i*n-n:i*n] = -B_sub
            A2[i*n:i*n+n,i*n:i*n+n] = B_sub

    # Factors
    X = (I + A1)
    Y = (I - A2)
    Z = (I - A1)
    W = (I + A2)
    # Inverse the required matrices
    X_inv = np.linalg.inv(X)
    W_inv = np.linalg.inv(W)
    # Final outputs
    P = np.matmul(X_inv,Y)
    Q = np.matmul(W_inv,Z)


    return P, Q

################################################################################
# Graphing function, handles graphing of schemes with exact solution
# Refactors z_vals 1-D input vector into 2-D array suitable for graphing
def graph(out,title,x_vals,y_vals,z_vals):
    # out    - title of output file
    # title  - graph title
    # x_vals - holds x axis values
    # y_vals - holds y axis values
    # z_vals - holds z values

    # X range limits
    x_l = min(x_vals)
    x_h = max(x_vals)
    y_l = min(y_vals)
    y_h = max(y_vals)

    # Refactor to 2D
    z_mesh = np.reshape(z_vals,(len(x_vals),len(y_vals)))

    # Plotting statements
    plt.contourf(x_vals,y_vals,z_mesh)
    plt.xlim((x_l,x_h))
    plt.ylim((y_l,y_h))
    plt.colorbar()
    plt.ylabel('Y')
    plt.xlabel('X')
    plt.title(title)
    plt.savefig(out)
    plt.clf()

################################################################################
# Main
################################################################################

print("Hello. Making Waves")

# Common variables
h = 1.0/30.0
lmda = 0.4
k = h*lmda

print(h,k,lmda)

# Initialize values
x_i = -1.0 # minimum x
x_f = 1.0 # maximum x
y_i = -1.0 # minimum y
y_f = 1.0 # maximum y
t_i = 0.0 # minimum time
t_f = 1.0 # maximum time

t_c = t_i

x_vals = []
y_vals = []
up_vals = []

# fill initial conditions - This method assumes a square domain and dx = dy
r_c = x_i
while r_c <= x_f:
    x_vals.append(r_c)
    y_vals.append(r_c)
    r_c += h

for valy in y_vals:
    for valx in x_vals:
        up_vals.append(u_exact(valx,valy))

N = len(x_vals) # Assumes num of grid points in x == grid points in y

# Fill matrices
P, Q = mat_fill(h,k,N)

r = 0
t_c = t_i
output = "grph_t"+str(r).zfill(3)+".jpg".format(r)
title = "ADI Solved System \n t = %.3f"% t_c
graph(output,title,x_vals,y_vals,up_vals)

# Update system
print("Beginning Time Evolution")
t_c += k
r += 1
while t_c <= t_f + k:

    u_inter = np.matmul(Q,up_vals)

    # uphold boundary conditions
    for j, valy in enumerate(y_vals):
        for i, valx in enumerate(x_vals):
            if valy == y_i:
                u_inter[N*j+i] = u_exact(valx-(t_c+k/2),valy-(t_c+k/2))
            if valx == x_i:
                u_inter[N*j+i] = u_exact(valx-(t_c+k/2),valy-(t_c+k/2))
            if valy == y_f:
                u_inter[N*j+i] = up_vals[N*(i-1)+i]
            if valx == x_f:
                u_inter[N*j+i] = up_vals[N*i+i-1]

    u_n = np.matmul(P,u_inter)


    # format for output graphs
    output = "grph_t"+str(r).zfill(3)+".jpg".format(r)
    title = "ADI Solved System \n t = %.3f"% t_c
    #graph(output,title,x_vals,y_vals,u_n)

    up_vals = u_n[:]

    # exact solution with plotting
    u_ex = []
    for valy in y_vals:
        for valx in x_vals:
                loc_val = u_exact(valx-t_c,valy-2*t_c)
                u_ex.append(loc_val)

    ex_title = "exact_t"+str(r).zfill(3)+".jpg".format(r)
    titl_e = "Exact Solution \n t = %.3f"% t_c
    #graph(ex_title,titl_e,x_vals,y_vals,u_ex)

    u_ex.clear()
    t_c += k
    r += 1

# Creates a gif of the outputs. comment out if unwanted.
#os.system("ffmpeg -y -i 'grph_t%03d.jpg' num_sln.gif")
#os.system("ffmpeg -y -i 'exact_t%03d.jpg' exact_sln.gif")
#os.system("rm -f *.jpg")
################################################################################
