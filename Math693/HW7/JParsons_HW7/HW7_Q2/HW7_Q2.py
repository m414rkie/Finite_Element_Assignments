#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Program written for Homework 7, MATH 693B. Problem outline from
# Strikwerda 10.4.1

# Program solves the heat equation for x in [-1,1] and t in [0,1]. Periodic
# boundaries will be used in conjunction with the Forward-Time Central-Space
# explicit scheme. This program will demonstrate second order accuracy using
# h = [1/10,1/20,1/40,1/80] as well as two initial conditions

# Author: Jon Parsons
# Date: 5-5-20

################################################################################
#Functions
################################################################################
# Function containing the exact solution
def u_ex_a(x,t):
    # x - point of interest
    # t - current time

    l = 50
    pi = np.pi

    i = 0
    sum = 0.0
    while i <= l:
        kk = 2*i+1
        sum += np.exp(-t*(kk*pi)**2)*(((-1)**i)/kk)*np.cos(kk*pi*x)
        i += 1


    sln = 0.5 + (2/pi)*sum

    return sln

################################################################################
# Function containing the exact solution
def u_ex_b(x,t):
    # x - point of interest
    # t - current time

    pi = np.pi

    sln = np.cos(pi*x)*np.exp(-pi*pi*t)

    return sln

################################################################################
# Function contains the first set of initial conditions
def ini_a(x):
    # x - point of interest

    if abs(x) < 0.5:
        sln = 1
    elif abs(x) == 0.5:
        sln = 0.5
    else:
        sln = 0

    return sln

################################################################################
# Function contains the second set of initial conditions
def ini_b(x):
    # x - point of interest

    sln = np.cos(np.pi*x)

    return sln

################################################################################
# Function containing the FTCS scheme
def FTCS(u_nm,u_nmx,u_nmp,mu):
    # u_nm  - Value at local point, current time
    # u_nmx - Value at m+1, current time
    # u_nmp - Value at m-1, current time
    # mu    - k/h^2


    u_n1 = u_nm + mu*(u_nmx - 2*u_nm + u_nmp)

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
# Function contains code to find the order of accuracy.
def acc_ord(v):
    # v - vector containing the error from the norms

    sum = 0
    N = 0
    while N < len(v)-1:
        ord_c = np.log2(v[N]/v[N+1])
        sum += ord_c
        print(ord_c)
        N += 1

    print(N)
    ord = sum/(N)

    return ord

################################################################################
# Graphing function, handles graphing of schemes with exact solution
def graph(out,title,x1_vals,y_vals,y2_vals):
    # out     - title of output file
    # title   - graph title
    # x_vals  - holds x axis values
    # y_vals  - holds B1 calculated values
    # y2_vals - holds exact values

    # X range limits
    x_l = min(x_vals)
    x_h = max(x_vals)

    # y range limits.
    y_h = 1.2
    y_l = -1.25


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

print("Hello. Converging waves.")

# Common variables
h_vals = [1/10,1/20,1/40,1/80]
mu = 0.4

# Initialize values
x_i = -1.0 # minimum x
x_f = 1.0 # maximum x
t_i = 0.0 # minimum time
t_f = 1.0 # maximum time

I1_l2_norms = np.zeros(len(h_vals))
I2_l2_norms = np.zeros(len(h_vals))

q = 0
for h in h_vals:
    q += 1
    k = h*h*mu
    N = round((x_f-x_i)/h) -1
    print("Solving for h = ",h)

    # Initialize vectors for previous values
    x_vals = np.zeros(N+1)
    I1_p   = np.zeros(N+1)
    I2_p   = np.zeros(N+1)
    I1_ex  = np.zeros(N+1)
    I2_ex  = np.zeros(N+1)
    # Initialize vectors for the new values
    I1_n   = np.zeros(N+1)
    I2_n   = np.zeros(N+1)

    # fill initial conditions
    x_c = x_i
    i = 0
    while i <= N:
        x_vals[i] = x_c
        I1_ex[i] = u_ex_a(x_c,t_i)
        I2_ex[i] = u_ex_b(x_c,t_i)
        I1_p[i]  = u_ex_a(x_c,t_i)
        I2_p[i]  = u_ex_b(x_c,t_i)

        x_c += h
        i += 1

    # Output graph of initial conditions
    t_c = t_i + k
    r = 0
    output = "grph_t"+str(r).zfill(3)+"h"+str(q)+"i1.jpg".format(r,q)
    title = "Waves With Boundaries \n t = %.3f"% t_c
    graph(output,title,x_vals,I1_p,I1_ex)
    output = "grph_t"+str(r).zfill(3)+"h"+str(q)+"i2.jpg".format(r,q)
    title = "Waves With Boundaries \n t = %.3f"% t_c
    graph(output,title,x_vals,I2_p,I2_ex)

    # Update system
    t_c += k
    r += 1
    while t_c <= t_f:
        # Interior points
        j = 0
        while j <= N:
            if j == 0:
                # Left Boundaries
                I1_n[j]  = FTCS(I1_p[j],I1_p[j+1],I1_p[N],mu)
                I2_n[j]  = FTCS(I2_p[j],I2_p[j+1],I2_p[N],mu)
                I1_ex[j] = u_ex_a(x_i,t_c)
                I2_ex[j] = u_ex_b(x_i,t_c)

            elif j == N:
                # Left Boundaries
                I1_n[j]  = FTCS(I1_p[j],I1_p[0],I1_p[j-1],mu)
                I2_n[j]  = FTCS(I2_p[j],I2_p[0],I2_p[j-1],mu)
                I1_ex[j] = u_ex_a(x_i,t_c)
                I2_ex[j] = u_ex_b(x_i,t_c)

            else:
                I1_n[j]  = FTCS(I1_p[j],I1_p[j+1],I1_p[j-1],mu)
                I2_n[j]  = FTCS(I2_p[j],I2_p[j+1],I2_p[j-1],mu)
                I1_ex[j] = u_ex_a(x_vals[j],t_c)
                I2_ex[j] = u_ex_b(x_vals[j],t_c)

            j += 1

        # Graph
        output = "grph_t"+str(r).zfill(3)+"h"+str(q)+"i1.jpg".format(r,q)
        title = "Waves With Boundaries \n t = %.3f"% t_c
        graph(output,title,x_vals,I1_p,I1_ex)

        output = "grph_t"+str(r).zfill(3)+"h"+str(q)+"i2.jpg".format(r,q)
        title = "Waves With Boundaries \n t = %.3f"% t_c
        graph(output,title,x_vals,I2_p,I2_ex)

        # Update lists
        for i in range(len(I1_n)):
            I1_p[i] = I1_n[i]
            I2_p[i] = I2_n[i]

        t_c += k
        r += 1

    I1_l2_norms[q-1] = L2_norm(I1_p,I1_ex,h)
    I2_l2_norms[q-1] = L2_norm(I2_p,I2_ex,h)

print("Determining Order of Accuracy")

print(I1_l2_norms)
print(I2_l2_norms)


I1_ord = acc_ord(I1_l2_norms)
I2_ord = acc_ord(I2_l2_norms)

print("Order for I1: ", I1_ord)
print("Order for I2: ", I2_ord)

# Creates a gif of the outputs. comment out if unwanted.
os.system("ffmpeg -y -i 'grph_t%03dh2i1.jpg' coupled_eqns1.gif")
os.system("ffmpeg -y -i 'grph_t%03dh2i2.jpg' coupled_eqns2.gif")

#os.system("rm -f *.jpg")
os.system("rm -f *h1*.jpg")
os.system("rm -f *h2*.jpg")
os.system("rm -f *h3*.jpg")



################################################################################
