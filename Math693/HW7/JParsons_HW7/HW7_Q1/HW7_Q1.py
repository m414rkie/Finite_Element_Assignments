#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Program written for Homework 7, MATH 693B. Problem outline from
# Strikwerda 10.3.6

# Program solves the first order wave equation on x in [-1.1], t in [0,0.96]
# We will use the Lax-Wendroff scheme with 4 values of h, (1/10,1/20,1/40,1/80)
# Three initial conditions will be used and compared via rate of convergance.
# Boundaries will be periodic

# Author: Jon Parsons
# Date: 5-5-20

################################################################################
#Functions
################################################################################
# Function containing the first initial conditions, as well as exact solution
def u_ex_a(x,t):
    # x - point of interest
    # t - current time

    if abs(x-t) < 0.5:
        sln = 1
    elif abs(x-t) == 0.5:
        sln = 0.5
    else:
        sln = 0

    return sln

################################################################################
# Function containing the second initial conditions, as well as exact solution
def u_ex_b(x,t):
    # x - point of interest
    # t - current time

    sln = np.cos(np.pi*(x-t))

    return sln

################################################################################
# Function containing the third initial conditions, as well as exact solution
def u_ex_c(x,t):
    # x - point of interest
    # t - current time

    if abs(x-t) <= 0.5:
        sln = (np.cos(np.pi*(x-t)))**2
    else:
        sln = 0

    return sln

################################################################################
# Function containing the LW scheme
def LW(u_nm,u_nmx,u_nmp,k,h):
    # u_nm  - Value at local point, current time
    # u_nmx - Value at m+1, current time
    # u_nmp - Value at m-1, current time
    # k     - Stepsize in time
    # h     - Stepsize in space

    lmda = k/h

    u_n1 = u_nm - 0.5*lmda*(u_nmx-u_nmp) + 0.5*lmda*lmda*(u_nmx-2*u_nm+u_nmp)

    return u_n1

################################################################################
# Function contains code to find the L2 norm
def L2_norm(v,e,h):
    # v - Vector containing the numerical solution
    # e - Vector containing the exact solution
    # h - stepsize

    sum = 0
    for i in range(len(v)):
        sum += (e[i] - v[i])**2

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
        ord_c = np.log2(v[N]/v[N+1])
        sum += ord_c
        print(ord_c)
        N += 1

    ord = sum/(N)

    return ord

################################################################################
# Graphing function, handles graphing of schemes with exact solution
def graph(out,title,x1_vals,y_vals,y2_vals):
    # out     - title of output file
    # title   - graph title
    # x_vals  - holds x axis values
    # y_vals  - holds calculated values
    # y2_vals - holds exact values

    # X range limits
    x_l = min(x_vals)
    x_h = max(x_vals)

    # y range limits.
    y_h = 2.0
    y_l = -2.0


    plt.plot(x_vals,y_vals)
    plt.plot(x_vals,y2_vals)
    plt.legend(["Numerical","Exact"],\
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
lmda = 0.8

# Initialize values
x_i = -1.0 # minimum x
x_f = 1.0 # maximum x
t_i = 0.0 # minimum time
t_f = 0.96 # maximum time

I1_l2_norms = np.zeros(len(h_vals))
I2_l2_norms = np.zeros(len(h_vals))
I3_l2_norms = np.zeros(len(h_vals))

q = 0
for h in h_vals:
    q += 1
    k = h*lmda
    N = round((x_f-x_i)/h)
    print("Solving for h = ",h)
    print(lmda,k,h)

    # Initialize vectors for holding previous values (And the exact vector)
    x_vals = np.zeros(N+1)
    I1_p   = np.zeros(N+1)
    I2_p   = np.zeros(N+1)
    I3_p   = np.zeros(N+1)
    I1_Ex  = np.zeros(N+1)
    I2_Ex  = np.zeros(N+1)
    I3_Ex  = np.zeros(N+1)
    # Initialize vectors for the next values
    I1_n   = np.zeros(N+1)
    I2_n   = np.zeros(N+1)
    I3_n   = np.zeros(N+1)

    # fill initial conditions
    x_c = x_i
    i = 0
    while i <= N:
        x_vals[i] = x_c
        I1_p[i]   = u_ex_a(x_c,t_i)
        I2_p[i]   = u_ex_b(x_c,t_i)
        I3_p[i]   = u_ex_c(x_c,t_i)
        I1_Ex[i]  = u_ex_a(x_c,t_i)
        I2_Ex[i]  = u_ex_b(x_c,t_i)
        I3_Ex[i]  = u_ex_c(x_c,t_i)
        x_c += h
        i += 1

    # Output graph of initial conditions
    t_c = t_i + k
    r = 0
    output = "grph_t"+str(r).zfill(3)+"h"+str(q)+"i1.jpg".format(r,q)
    title = "Waves With Periodic Boundaries \n t = %.3f"% t_c
    graph(output,title,x_vals,I1_p,I1_Ex)

    output = "grph_t"+str(r).zfill(3)+"h"+str(q)+"i2.jpg".format(r,q)
    title = "Waves With Periodic Boundaries \n t = %.3f"% t_c
    graph(output,title,x_vals,I2_p,I2_Ex)

    output = "grph_t"+str(r).zfill(3)+"h"+str(q)+"i3.jpg".format(r,q)
    title = "Waves With Periodic Boundaries \n t = %.3f"% t_c
    graph(output,title,x_vals,I3_p,I3_Ex)

    # Update system
    t_c += k
    r += 1
    while t_c <= t_f:

        j = 0
        while j <= N:
            if j == 0:
                # Periodic Boundaries
                I1_n[0]  = LW(I1_p[0],I1_p[1],I1_p[N],k,h)
                I2_n[0]  = LW(I2_p[0],I2_p[1],I2_p[N],k,h)
                I3_n[0]  = LW(I3_p[0],I3_p[1],I3_p[N],k,h)
                I1_Ex[j] = u_ex_a(x_i,t_c)
                I2_Ex[j] = u_ex_b(x_i,t_c)
                I3_Ex[j] = u_ex_c(x_i,t_c)
            elif j == N:
                I1_n[N]  = LW(I1_p[N],I1_p[0],I1_p[N-1],k,h)
                I2_n[N]  = LW(I2_p[N],I2_p[0],I2_p[N-1],k,h)
                I3_n[N]  = LW(I3_p[N],I3_p[0],I3_p[N-1],k,h)
                I1_Ex[j] = u_ex_a(x_f,t_c)
                I2_Ex[j] = u_ex_b(x_f,t_c)
                I3_Ex[j] = u_ex_c(x_f,t_c)
            else:
                I1_n[j]  = LW(I1_p[j],I1_p[j+1],I1_p[j-1],k,h)
                I2_n[j]  = LW(I2_p[j],I2_p[j+1],I2_p[j-1],k,h)
                I3_n[j]  = LW(I3_p[j],I3_p[j+1],I3_p[j-1],k,h)
                I1_Ex[j] = u_ex_a(x_vals[j],t_c)
                I2_Ex[j] = u_ex_b(x_vals[j],t_c)
                I3_Ex[j] = u_ex_c(x_vals[j],t_c)
            j += 1

        # Graph
        output = "grph_t"+str(r).zfill(3)+"h"+str(q)+"i1.jpg".format(r,q)
        title = "Waves With Periodic Boundaries \n t = %.3f"% t_c
        graph(output,title,x_vals,I1_n,I1_Ex)

        output = "grph_t"+str(r).zfill(3)+"h"+str(q)+"i2.jpg".format(r,q)
        title = "Waves With Periodic Boundaries \n t = %.3f"% t_c
        graph(output,title,x_vals,I2_n,I2_Ex)

        output = "grph_t"+str(r).zfill(3)+"h"+str(q)+"i3.jpg".format(r,q)
        title = "Waves With Periodic Boundaries \n t = %.3f"% t_c
        graph(output,title,x_vals,I3_n,I3_Ex)

        # Update lists
        for i in range(len(I1_p)):
            I1_p[i] = I1_n[i]
            I2_p[i] = I2_n[i]
            I3_p[i] = I3_n[i]

        t_c += k
        r += 1

    I1_l2_norms[q-1] = L2_norm(I1_p,I1_Ex,h)
    I2_l2_norms[q-1] = L2_norm(I2_p,I2_Ex,h)
    I3_l2_norms[q-1] = L2_norm(I3_p,I3_Ex,h)

print("Determining Order of Accuracy")

print(I1_l2_norms)
print(I2_l2_norms)
print(I3_l2_norms)

I1_ord = acc_ord(I1_l2_norms)
I2_ord = acc_ord(I2_l2_norms)
I3_ord = acc_ord(I3_l2_norms)

print("Order for I1: ", I1_ord)
print("Order for I2: ", I2_ord)
print("Order for I3: ", I3_ord)

# Creates a gif of the outputs. comment out if unwanted.
os.system("ffmpeg -y -i 'grph_t%03dh4i1.jpg' coupled_eqns1.gif")
os.system("ffmpeg -y -i 'grph_t%03dh4i2.jpg' coupled_eqns2.gif")
os.system("ffmpeg -y -i 'grph_t%03dh4i3.jpg' coupled_eqns3.gif")

#os.system("rm -f *.jpg")
#os.system("rm -f *h1*.jpg")
#os.system("rm -f *h2*.jpg")
#os.system("rm -f *h3*.jpg")


################################################################################
