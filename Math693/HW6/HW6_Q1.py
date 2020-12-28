#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Program written for Homework 6, MATH 693B. Problem outline from
# Strikwerda 8.3.1

# Program solves the wave equation u_tt = u_xx with a variety of values for h
# and lambda held at one. We will use the second order Central-Time
# Central-Space (CTCS) scheme with three sets of boundary conditions, two first
# order accurate, the other second order accurate. These will be denoted as
# B1 and B2 and B3 respectively. In order to preserve the properties of the
# boundary conditions without degradation, the timestep '-1' will be the exact
# solution

# Author: Jon Parsons
# Date: 4-15-20

################################################################################
#Functions
################################################################################
# Function containing the initial condition of ODE u, as well as exact solution
def u_ex(x,t):
    # x - point of interest
    # t - current time

    sln = np.cos(x + t) + np.cos(x - t)

    return sln

################################################################################
# Function containing the CTCS scheme
def CTCS(u_nm,u_nmx,u_nmp,u_npm,k,h):
    # u_nm  - Value at local point, current time
    # u_npm - Value at local point, previous time
    # u_nmx - Value at m+1, current time
    # u_nmp - Value at m-1, current time
    # k     - Stepsize in time
    # h     - Stepsize in space

    C = (k/h)**2

    u_n1 = 2*u_nm - u_npm + C*(u_nmx - 2*u_nm + u_nmp)

    return u_n1

################################################################################
# Function containing the B1 conditions, a second order accurate scheme
def B1(u_nxmx,u_nxmxx):
    # u_nxm   - Value at point m+1, next time
    # u_nxmxx - Value at point m+2, next time

    u_n1 = (4.0*u_nxmx-u_nxmxx)/3.0

    return u_n1

################################################################################
# Function containing the B2 conditions, a second order accurate scheme
def B2(u_nm,u_npm,u_nmx,k,h):
    # u_nm  - Value at local point
    # u_npm - Value at local point, previous time
    # u_nmx - Value at point m+1, current time
    # k     - Stepsize in time
    # h     - Stepsize in space

    C = (k/h)**2
    u_n1 = 2.0*u_nm - u_npm - 2*C*(u_nm-u_nmx)

    return u_n1

################################################################################
# Function containing the B3 conditions, a first order accurate scheme
def B3(u_nmx):
    # u_nmx - Value at point m+1, current time

    u_n1 = u_nmx

    return u_n1
################################################################################
# Function contains code to find the L2 norm
def L2_norm(v,e,h):
    # v - Vector containing the numerical solution
    # e - Vector containing the exact solution
    # h - stepsize
    sum = 0
    for i in range(len(v)):
        sum += (e[i] - v[i])*(e[i] - v[i])

    l2 = np.sqrt(sum*h)

    return l2

################################################################################
# Function contains code to find the order of accuracy. Arguements must be of
# equal length
def acc_ord(v):
    # v - vector containing the error from the norms

    sum = 0
    N = 0
    while N < len(v)-1:
        ord_c = np.log(v[N]/v[N+1])/np.log(2)
        sum += ord_c
        print(ord_c)
        N += 1

    print(N)
    ord = sum/(N)

    return ord

################################################################################
# Graphing function, handles graphing of schemes with exact solution
def graph(out,title,x1_vals,y_vals,y2_vals,y3_vals,y4_vals):
    # out     - title of output file
    # title   - graph title
    # x_vals  - holds x axis values
    # y_vals  - holds B1 calculated values
    # y2_vals - holds B2 calculated values
    # y3_vals - holds B3 calculated values
    # y4_vals - holds exact values

    # X range limits
    x_l = min(x_vals)
    x_h = max(x_vals)

    # y range limits.
    y_h = 2.1
    y_l = 0.0


    plt.plot(x_vals,y_vals)
    plt.plot(x_vals,y2_vals)
    plt.plot(x_vals,y3_vals)
    plt.plot(x_vals,y4_vals)
    plt.legend(["2nd Order (B1)","2nd Order (B2)","1st Order (B3)","Exact"],\
                loc='upper right')
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

print("Hello. Expanding waves.")

# Common variables
h_vals = [1/10,1/20,1/40,1/80]
lmda = 1.0

# Initialize values
x_i = 0.0 # minimum x
x_f = 1.0 # maximum x
t_i = 0.0 # minimum time
t_f = 1.0 # maximum time

B1_l2_norms = []
B2_l2_norms = []
B3_l2_norms = []

q = 0
for h in h_vals:
    q += 1
    k = h*lmda
    N = round((x_f-x_i)/h)
    print("Solving for h = ",h)

    # Initialize vectors for holding previous values (And the exact vector)
    x_vals = np.zeros(N+1)
    B1_p   = np.zeros(N+1)
    B2_p   = np.zeros(N+1)
    B3_p   = np.zeros(N+1)
    Ex     = np.zeros(N+1)
    # Initialize vectors for current values
    B1_c   = np.zeros(N+1)
    B2_c   = np.zeros(N+1)
    B3_c   = np.zeros(N+1)
    # Initialize vectors for the next values
    B1_n   = np.zeros(N+1)
    B2_n   = np.zeros(N+1)
    B3_n   = np.zeros(N+1)

    # fill initial conditions
    x_c = x_i
    i = 0
    while i <= N:
        x_vals[i] = x_c
        Ex[i] = u_ex(x_c,t_i)
        B1_p[i] = u_ex(x_c,t_i)
        B2_p[i] = u_ex(x_c,t_i)
        B3_p[i] = u_ex(x_c,t_i)
        B1_c[i] = u_ex(x_c,t_i+k)
        B2_c[i] = u_ex(x_c,t_i+k)
        B3_c[i] = u_ex(x_c,t_i+k)
        x_c += h
        i += 1

    # Output graph of initial conditions
    t_c = t_i + k
    r = 0
    output = "grph_t"+str(r).zfill(3)+"h"+str(q)+".jpg".format(r,q)
    title = "Waves With Boundaries \n t = %.3f"% t_c
    graph(output,title,x_vals,B1_c,B2_c,B3_c,Ex)

    # Update system
    t_c += k
    r += 1
    while t_c <= t_f:
        # Interior points
        j = N
        while j >= 0:
            if j == N:
                # Right boundaries
                B1_n[j] = u_ex(x_f,t_c)
                B2_n[j] = u_ex(x_f,t_c)
                B3_n[j] = u_ex(x_f,t_c)
                Ex[j]   = u_ex(x_f,t_c)
            elif j == 0:
                # Left Boundaries
                B1_n[j] = B1(B1_n[j+1],B1_n[j+2])
                B2_n[j] = B2(B2_c[j],B2_p[j],B2_c[j+1],k,h)
                B3_n[j] = B3(B3_n[j+1])
                Ex[j]   = u_ex(x_i,t_c)
            else:
                B1_n[j] = CTCS(B1_c[j],B1_c[j+1],B1_c[j-1],B1_p[j],k,h)
                B2_n[j] = CTCS(B2_c[j],B2_c[j+1],B2_c[j-1],B2_p[j],k,h)
                B3_n[j] = CTCS(B3_c[j],B3_c[j+1],B3_c[j-1],B3_p[j],k,h)
                Ex[j]   = u_ex(x_vals[j],t_c)
            j -= 1

        # Graph
        output = "grph_t"+str(r).zfill(3)+"h"+str(q)+".jpg".format(r,q)
        title = "Waves With Boundaries \n t = %.3f"% t_c
        graph(output,title,x_vals,B1_n,B2_n,B3_n,Ex)
        # Update lists
        for i in range(len(B1_c)):
            B1_p[i] = B1_c[i]
            B2_p[i] = B2_c[i]
            B3_p[i] = B3_c[i]
            B1_c[i] = B1_n[i]
            B2_c[i] = B2_n[i]
            B3_c[i] = B3_n[i]

        t_c += k
        r += 1

    B1_l2_norms.append(L2_norm(B1_c,Ex,h))
    B2_l2_norms.append(L2_norm(B2_c,Ex,h))
    B3_l2_norms.append(L2_norm(B3_c,Ex,h))

print("Determining Order of Accuracy")

print(B1_l2_norms)
print(B2_l2_norms)
print(B3_l2_norms)

B1_ord = acc_ord(B1_l2_norms)
B2_ord = acc_ord(B2_l2_norms)
B3_ord = acc_ord(B3_l2_norms)

print("Order for B1: ", B1_ord)
print("Order for B2: ", B2_ord)
print("Order for B3: ", B3_ord)

# Creates a gif of the outputs. comment out if unwanted.
#os.system("ffmpeg -y -i 'grph_t%03dh1.jpg' coupled_eqns1.gif")
#os.system("ffmpeg -y -i 'grph_t%03dh2.jpg' coupled_eqns2.gif")
#os.system("ffmpeg -y -i 'grph_t%03dh3.jpg' coupled_eqns3.gif")

#os.system("rm -f *.jpg")



################################################################################
