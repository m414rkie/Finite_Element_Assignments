#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Program written for Homework 4, MATH 693B. Problem outline from
# Strikwerda 3.2.4

# Program solves a driven ODE using the box method. Order of accuracy is checked
# using both the L2 and max norm.

# Author: Jon Parsons
# Date: 3-29-20

################################################################################
#Functions
################################################################################
# Function containing the initial conditions
def ini_fill(x):
    # x - point of interest

    sln = np.sin(x)

    return sln

################################################################################
# Function containing the driving function of the system
def drvr(x,t):
    # x - point of interest
    # t - time of the system

    sln = np.sin(x-t)

    return sln

################################################################################
# Function containing the left boundary condition
def l_bnd(t):
    # t - time of the system

    sln = -(1.0 + t)*np.sin(t)

    return sln

################################################################################
# Function containing the update eqn of u
def u_updt(f2_p1,f2_nn,f_p1,f_nn,u2_nn,u_nn,u_p1,h,k):
    # f2_p1 - driving function, current time at point of interest
    # f2_nn - driving function, current time left of point of interest
    # f_p1  - driving function, previous time at point of interest
    # f_nn  - driving function, previous time left of point of interest
    # u2_nn - wave, current time left of point of interest
    # u_nn  - wave, previous time left of point of interest
    # u_p1  - wave, previous time at point of interest
    # h     - dimension size
    # k     - timestep size

    # dividing constant
    lmda = k/h

    # driving portion
    drvr = (k/(2.0*(1+lmda)))*(f2_p1+f2_nn+f_p1+f_nn)

    # combined
    u_n1 = drvr + u_nn + ((1-lmda)/(1+lmda))*(u_p1-u2_nn)

    return u_n1

################################################################################
# Function containing code for L2 norm
def L2_norm(arr,dx):
    # arr - vector to find the norm for

    # initialize values
    sum = 0.0
    for val in arr:
        sum += val*val

    norm = np.sqrt(sum*dx)
    return norm

################################################################################
# Function containing code for the max norm
def Max_norm(arr):
    # arr - vector to find the norm for

    # Must be done this way, abs() cannot operate on vector
    norm = 0.0
    for val in arr:
        if (abs(val) > norm):
            norm = abs(val)

    return norm

################################################################################
# Function containing code for determining order of accuracy. If there are
#  values to give multiple determinations, an average is returned.
# Function assumes that h is halved at each iteration
def ord_acc(vals):
    # vals - vector containing norms, can be any length greater than 2

    if len(vals) < 2:
        print("Insufficient number of norms to calculate order of accuracy.")
        return 0

    sum = 0
    i = 1
    while i < len(vals) - 1:
        p_int = 0.0
        p_int = (vals[i-1] - vals[i])/(vals[i] - vals[i+1])
        print("p_int",np.log2(p_int))
        sum += np.log2(p_int)
        i +=1

    p = sum/(i-1)
    return p

################################################################################
# Graphing function, handles graphing of schemes with exact solution
def graph(out,title,x1_vals,y_vals,y2_vals):
    # out    - title of output file
    # title  - graph title
    # x_vals - holds x axis values
    # y_vals - holds resultant values
    # y2_vals- holds driving values

    # X range limits
    x_l = min(x_vals)
    x_h = max(x_vals)
    # y range limits.
    y_l = -5.0
    y_h = 5.0

    plt.plot(x_vals,y_vals)
    plt.plot(x_vals,y2_vals)
    plt.legend(["Wave","Driver"], loc='upper right')
    plt.xlim((x_l,x_h))
    plt.ylim((y_l,y_h))
    plt.ylabel('Wave Magnitude')
    plt.xlabel('X')
    plt.title(title)
    plt.savefig(out)
    plt.clf()

################################################################################
# Main
################################################################################

print("Hello. Driving Forward")

# Common variables
h_vals = [1.0/10.0,1.0/20.0,1.0/40.0,1.0/80.0]
lmda = 1.2

# Initialize values
x_i = 0.0 # minimum x
x_f = 1.0 # maximum x
t_i = 0.0 # minimum time
t_f = 1.2 # maximum time, extended due to gif maker

# Initialize arrays for holding norms
l2norm_vals = []
mxnorm_vals = []

# Outer loop iterates through h values
for h in h_vals:
    # Determine step size in time
    k = lmda*h
    # Initialize arrays, u is wave, d is driver. p is previous, n is current
    x_vals  = []
    u_vals = []
    d_vals = []
    u2_vals = []
    d2_vals = []
    # Start at left side of system
    x_c = x_i
    t_c = t_i

    # fill initial conditions
    while x_c <= x_f:
        x_vals.append(x_c)
        u_vals.append(ini_fill(x_c))
        d_vals.append(drvr(x_c,t_c))
        x_c += h

    # Graph initial conditions
    r = 0
    output = "grph_t"+str(r).zfill(3)+".jpg".format(r)
    title = ("Coupled System \n t = {:.3f} h = {:.2f}".format(t_c,h))
    graph(output,title,x_vals,u_vals,d_vals)

    N = len(u_vals)

    # Update system
    r += 1
    while t_c <= t_f + k:
        # Update time
        t_c += k
        # Fill driving values
        for x in x_vals:
            d2_vals.append(drvr(x,t_c))

        # Boundary conditions for u, left
        u2_vals.append(l_bnd(t_c))

        for i in range(1,N):
            u2_vals.append(u_updt(d2_vals[i],d2_vals[i-1],d_vals[i],d_vals[i-1], \
                         u2_vals[i-1],u_vals[i-1],u_vals[i],h,k))

        # format for output graphs
        output = "grph_t"+str(r).zfill(3)+".jpg".format(r)
        title = ("Coupled System \n t = {:.3f} h = {:.2f}".format(t_c,h))
        graph(output,title,x_vals,u2_vals,d2_vals)

        # Update and reset
        u_vals = u2_vals[:]
        d_vals = d2_vals[:]
        u2_vals.clear()
        d2_vals.clear()

        r += 1

    l2norm_vals.append(L2_norm(u_vals,h))
    mxnorm_vals.append(Max_norm(u_vals))

    # Creates a gif of the outputs. comment out if unwanted.
    #os.system("ffmpeg -y -i 'grph_t%03d.jpg' coupled_eqns_h{}.gif".format(h))
    #os.system("rm -f *.jpg")

print("Order of accuracy from L2 norm:")
l2_ord = ord_acc(l2norm_vals)
print("L2 norm values: ",l2norm_vals)
print(l2_ord)
print("Order of accuracy from Max norm:")
mx_ord = ord_acc(mxnorm_vals)
print("Max norm values: ",mxnorm_vals)
print(mx_ord)

print

################################################################################
