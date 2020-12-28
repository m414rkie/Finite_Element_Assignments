#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Program written for final project, MATH 693B.

# Program contains functions resulting from the von Neumann stability analysis
# of two explicit and two implicit finite elements methods. Program will output
# graphs of these on the complex plane
# Explicit: CTCS ; 2-step scheme
# Implicit: FTFS ; Crank-Nicolson

# Author: Jon Parsons
# Date: 4-20-20

################################################################################
#Functions
################################################################################
# Function containing the resulting equation of the von Neumann analysis of the
# FTFS scheme
def FTCS_stab(a,mu,om):
    # a    - wavespeed
    # lmda - k/h^2
    # om   - angle

    #sln = 1 + 2*a*mu*(np.cos(om) - 1)
    sln = 1 - 4*a*mu*(np.sin(om/2)**2)

    return sln

################################################################################
# Function containing the resulting equation of the von Neumann analysis of the
# Crank-Nicolson scheme
def CN_stab(a,mu,om):
    # a    - wavespeed
    # lmda - k/h^2
    # om   - angle

    c = a*mu*(np.cos(om)-1)

    sln = (c+1)/(1-c)

    return sln

################################################################################
# Function containing the resulting equation of the von Neumann analysis of the
# CTCS scheme
def CTCS_stab(a,lmda,om):
    # a    - wavespeed
    # lmda - k/h
    # om   - angle

    b = 1 - a*a*lmda*lmda*(np.sin(om/2)**2)

    slnp = ((np.sqrt(b + 0j)) + 1j*a*lmda*np.sin(om/2))**2
    slnm = ((np.sqrt(b + 0j)) - 1j*a*lmda*np.sin(om/2))**2

    return slnp, slnm

################################################################################
# Function containing the resulting equation of the von Neumann analysis of the
# Second Order scheme
def Ord2_stab(a,lmda,om):
    # a    - wavespeed
    # lmda - k/h
    # om   - angle

    a = 1 + a*a*lmda*lmda*(np.sin(om/2)**2)
    sq_arg = 1 - a*a

    slnp = (1 + np.sqrt(sq_arg + 0j))/(a)
    slnm = (1 - np.sqrt(sq_arg + 0j))/(a)

    return slnp, slnm

################################################################################
# Graphing function, handles graphing of schemes with exact solution
def graph(out,title,r_vals,i_vals,ucx,ucy):
    # out     - title of output file
    # title   - graph title
    # r_vals  - holds real
    # i_vals  - holds the imaginary components
    # ucx     - holds x component of unit circle
    # ucy     - holds y component of unit circle

    # X range limits
    x_l = -3
    x_h = 3

    # y range limits. Will auto-adjust based on lowest value of the exact sln
    y_l = -3
    y_h = 3

    plt.plot(r_vals,i_vals,'r')
    plt.plot(ucx,ucy,'b')
    plt.legend(["Stability Results","Unit Circle"], loc='upper right')
    plt.xlim((x_l,x_h))
    plt.ylim((y_l,y_h))
    plt.ylabel('Im')
    plt.xlabel('Re')
    plt.title(title)
    plt.savefig(out)
    plt.clf()

################################################################################
# Graphing function, handles graphing. This version for quadratic
def graph_qd(out,title,r_vals,i_vals,r2_vals,i2_vals,ucx,ucy):
    # out     - title of output file
    # title   - graph title
    # r_vals  - holds real for positive
    # i_vals  - holds the imaginary components for positive
    # r_vals  - holds real for negative
    # i_vals  - holds the imaginary components for negatice
    # ucx     - holds x component of unit circle
    # ucy     - holds y component of unit circle

    # X range limits
    x_l = -3
    x_h = 3

    # y range limits. Will auto-adjust based on lowest value of the exact sln
    y_l = -3
    y_h = 3

    plt.plot(r_vals,i_vals,'r')
    plt.plot(r2_vals,i2_vals,'r')
    plt.plot(ucx,ucy,'b')
    plt.legend(["Stability Results","","Unit Circle"], loc='upper right')
    plt.xlim((x_l,x_h))
    plt.ylim((y_l,y_h))
    plt.ylabel('Im')
    plt.xlabel('Re')
    plt.title(title)
    plt.savefig(out)
    plt.clf()

################################################################################
# Main
################################################################################

print("Hello. Seeing Stability.")

a_vals = [2,1,0.5]
lmda_vals = [2,1,0.5]

n = 1000
d_angle = 10*np.pi/n
angle_vals = np.zeros(n+1)
results = np.zeros(n+1,dtype=complex)
results2 = np.zeros(n+1,dtype=complex)

# Get the unit circle
ucx = np.zeros(n)
ucy = np.zeros(n)
i = 0
while i < n:
    a = i*2*np.pi/n
    ucx[i] = np.cos(a)
    ucy[i] = np.sin(a)
    i += 1

i = 0
while i <= n:
    angle_vals[i] = i*d_angle
    i += 1

for val in a_vals:
    for lmda in lmda_vals:
        title = "FTCS \n Mu = {}, a = {}".format(val,lmda)

        output = "FTCS_mu{}_a{}.jpg".format(lmda,val)

        for l, ang in enumerate(angle_vals):
            results[l] = FTCS_stab(val,lmda,ang)

        graph(output,title,results.real,results.imag,ucx,ucy)


        title = "CN \n Mu = {}, a = {}".format(val,lmda)

        output = "CN_mu{}_a{}.jpg".format(lmda,val)

        for l, ang in enumerate(angle_vals):
            results[l] = CN_stab(val,lmda,ang)

        graph(output,title,results.real,results.imag,ucx,ucy)


        title = "CTCS \n Lambda = {}, a = {}".format(val,lmda)

        output = "CTCS_l{}_a{}.jpg".format(lmda,val)

        for l, ang in enumerate(angle_vals):
            results[l], results2[l] = CTCS_stab(val,lmda,ang)

        graph_qd(output,title,results.real,results.imag,results2.real,\
                                                     results2.imag,ucx,ucy)

        title = "2-Step \n Lambda = {}, a = {}".format(val,lmda)

        output = "step_l{}_a{}.jpg".format(lmda,val)

        for l, ang in enumerate(angle_vals):
            results[l], results2[l] = Ord2_stab(val,lmda,ang)

        graph_qd(output,title,results.real,results.imag,results2.real,\
                                                     results2.imag,ucx,ucy)
################################################################################
