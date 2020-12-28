#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Program written for Homework 4, MATH 693B. Problem outline from
# Strikwerda 1.3.2

# Program solves a system of coupled ODE's using the Lax-Friedrichs method.
# The system is solved on x in [-3,3] and t in [0,2] and using boundary
# conditions which are different for each ODE (u, and w).

# Author: Jon Parsons
# Date: 3-29-20

################################################################################
#Functions
################################################################################
# Function containing the initial condition of ODE u
def u_ini(x):
    # x - point of interest

    if 0.0 > 1.0 - abs(x):
        sln = 0.0
    else:
        sln = 1.0 - abs(x)

    return sln

################################################################################
# Function containing the initial condition of ODE w
def w_ini(x):
    # x - point of interest

    if 0.0 > 1.0 - 2.0*abs(x):
        sln = 0.0
    else:
        sln = 1.0 - 2.0*abs(x)

    return sln

################################################################################
# Function containing the update eqn of u
def u_updt(t,u_p1,u_m1,u_nn,w_p1,w_m1,lmbda,k):
    # t    - current time of the system
    # u_p1 - value of u immediately right of point of interest
    # u_m1 - value of u immediately left of point of interest
    # u_nn - Value of u at the point of interest
    # w_p1 - value of w immediately right of the point of interest
    # w_m1 - value of w immediately left of the point of interest
    # lmbda- value of lambda
    # k    - timestep size

    u_n1 = 0.5*(u_p1+u_m1)-(lmbda/6.0)*(t-2.0)*(u_p1-u_m1)-(lmbda/3.0)*(t+1.0)*(w_p1-w_m1)-(k/3.0)*u_nn

    return u_n1

################################################################################
# Function containing the update eqn of u
def w_updt(t,w_p1,w_m1,w_nn,u_p1,u_m1,lmbda,k):
    # t    - current time of the system
    # w_p1 - value of u immediately right of point of interest
    # w_m1 - value of u immediately left of point of interest
    # w_nn - Value of u at the point of interest
    # u_p1 - value of w immediately right of the point of interest
    # u_m1 - value of w immediately left of the point of interest
    # lmbda- value of lambda
    # k    - timestep size

    w_n1 = 0.5*(w_p1+w_m1)-(lmbda/6.0)*(t+1.0)*(u_p1-u_m1)-(lmbda/6.0)*(2.0*t+1.0)*(w_p1-w_m1)-(k/3.0)*w_nn

    return w_n1

################################################################################
# Graphing function, handles graphing of schemes with exact solution
def graph(out,title,x1_vals,y_vals,y2_vals):
    # out    - title of output file
    # title  - graph title
    # x_vals - holds x axis values
    # y_vals - holds u values
    # y2_vals- holds w values

    # X range limits
    x_l = min(x_vals)
    x_h = max(x_vals)
    # y range limits. Will auto-adjust based on lowest value. Max always 1
    if min(y_vals) < min(y2_vals):
        y_l = min(y_vals)
    else:
        y_l = min(y2_vals)

    y_h = 1.0

    plt.plot(x_vals,y_vals)
    plt.plot(x_vals,y2_vals)
    plt.legend(["U","W"], loc='upper right')
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

print("Hello. Solving Systems")

# Common variables
h = 1.0/20.0
lmda = 0.5
k = h*lmda

print(h,k,lmda)

# Initialize values
x_i = -3.0 # minimum x
x_f = 3.0 # maximum x
t_i = 0.0 # minimum time
t_f = 2.25 # maximum time, extended due to gif maker

t_c = t_i

x_vals = []
up_vals = []
wp_vals = []
un_vals = []
wn_vals = []

x_c = x_i

# fill initial conditions
while x_c <= x_f:
    x_vals.append(x_c)
    u_c = u_ini(x_c)
    w_c = w_ini(x_c)
    up_vals.append(u_c)
    wp_vals.append(w_c)
    x_c += h

N = len(up_vals)

# Update system
t_c = t_i
r = 0
while t_c <= t_f + k:

    # Boundary conditions for u, left
    un_vals = [0.0]

    for i in range(1,N-1):
        u_n = u_updt(t_c,up_vals[i+1],up_vals[i-1],up_vals[i],wp_vals[i+1],wp_vals[i-1],lmda,k)
        w_n = w_updt(t_c,wp_vals[i+1],wp_vals[i-1],wp_vals[i],up_vals[i+1],up_vals[i-1],lmda,k)
        un_vals.append(u_n)
        wn_vals.append(w_n)

    # Remaining boundary conditions
    un_vals.append(0.0)
    wn_vals.append(wn_vals[-1])
    wn_vals.insert(0,wn_vals[0])

    # format for output graphs
    output = "grph_t"+str(r).zfill(3)+".jpg".format(r)
    title = "Coupled System \n t = %.3f"% t_c
    graph(output,title,x_vals,un_vals,wn_vals)

    up_vals = un_vals[:]
    wp_vals = wn_vals[:]

    un_vals.clear()
    wn_vals.clear()

    t_c += k
    r += 1


# Creates a gif of the outputs. comment out if unwanted.
os.system("ffmpeg -y -i 'grph_t%03d.jpg' coupled_eqns.gif")
os.system("rm -f *.jpg")
################################################################################
