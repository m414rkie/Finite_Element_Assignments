#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math

# Program written for Homework 3, MATH 693B. Problem outline from
# Strikwerda 6.3.10

# Problem solves an initial boundary problem using the Crank-Nicolson method.
# Program will compare values of h (1/10,1/20,1/40),  at lambda = 1 as well as
# mu = 10. Error will be checked against the exact solution using the L2 norm
# and supremum norm

# Author: Jon Parsons
# Date: 3-8-20

################################################################################
#Functions
################################################################################
# Function containing the exact solution and the Boundary conditions. Solves
#    for each point individually
def ex_sln(x,l,t):
    # x - point of interest
    # l - used in summation
    # t - current time

    Pi = math.pi
    sum = 0
    for i in range(l+1):
        l_md = 2*i + 1
        sum += (-1**i)*(math.cos(Pi*l_md*x)/(Pi*(2*i+1)))*math.exp(-((Pi*l_md)**2)*t)

    sln = 0.5 + 2*sum
    return sln

################################################################################
# Function containing the initial conditions of the system
def ini_fill(dx,xi,n):
    # dx - step size in x
    # xi - Initial x value
    # xf - Final x value

    wv = []
    xv = []
    x = xi
    for i in range(round(n)):
        xv.append(x)
        if abs(x) < 0.5:
            wv.append(1.0)
        elif abs(x) == 0.5:
            wv.append(0.5)
        elif abs(x) > 0.5:
            wv.append(0.0)
        x += dx

    return xv, wv

################################################################################
# Function that fills the A and B arrays of the Crank-Nicolson method
def mat_fill(mu,n):
    # mu - k/h^2
    # n  - number of elements

    A_m = np.zeros((round(n),round(n)))
    B_m = np.zeros((round(n),round(n)))

    for i in range(round(n)):
        if i == 0:
            A_m[0,0] = 1
            B_m[0,0] = 1
        elif i == n-1:
            A_m[i,i] = 1
            B_m[i,i] = 1
        else:
            A_m[i,i-1] = -mu
            A_m[i,i]   = 1 + 2*mu
            A_m[i,i+1] = -mu
            B_m[i,i-1] = mu
            B_m[i,i]   = 1 - 2*mu
            B_m[i,i+1] = mu

    return A_m, B_m

################################################################################
# Function that calculates the L2 norm between 2 vectors
def l2_norm(a,b,h):
    # a - first vector
    # b - second vector
    # h - stepsize

    sum = 0.0
    for i in range(len(a)):
        sum += (h*h)*(b[i] - a[i])**2

    norm = np.sqrt(sum)

    return norm

################################################################################
# Function that calculates the max-norm between 2 vectors
def max_norm(a,b):
    # a - first vector
    # b - second vector

    mx = 0.0
    for i in range(len(a)):
        c = abs(a[i] - b[i])
        if c > mx:
            mx = c

    return mx

################################################################################
# Graphing function, handles graphing of schemes with exact solution
def graph(out,title,x_vals,y_vals,y2_vals):
    # out    - title of output file
    # title  - graph title
    # x_vals - holds x axis values
    # y_vals - holds y values of first set of data
    # y2_vals- holds y values of second set of data

    x_l = min(x_vals)
    x_h = max(x_vals)
    y_l = min(y_vals) - 0.1
    y_h = max(y_vals) + 0.1

    plt.plot(x_vals,y_vals)
    plt.plot(x_vals,y2_vals)
    plt.legend(["Numerical","Exact"], loc='upper right')
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

print("Hello. Making Waves.")
# Common variables
h_vals = [0.1,0.05,0.025] # 1/10, 1/20, 1/40
type_vals = [1,2] # 1 - mu is constant, 2 - lambda is constant

# Initialize values
l = 25 # Summation max for exact solution
x_i = -1.0 # minimum x
x_f = 1.0 # maximum x
t_i = 0.0 # minimum time
t_f = 0.5 # maximum time

# Loops
# Determine constant mu or lambda
for tp in type_vals:
    # Determine h value
    for h in h_vals:
        N  = (x_f - x_i)/h # number of gridpoints
        if tp == 1:
            print("Constant Mu")
            dt = 10*h*h # Mu version
        elif tp == 2:
            print("Constant lambda")
            dt = h # lambda version

        # determine mu, lambda
        lmbda = dt/h
        Mu = dt/(h*h)

        # initial fill
        x_vals, wv_i = ini_fill(h,x_i,N)
        A, B = mat_fill(Mu,N) # To change wavespeed Mu -> D*Mu, D is speed
        A_inv = np.linalg.inv(A)
        AB = np.matmul(A_inv,B)

        # iterate through time
        t_c = t_i
        while t_c <= t_f:
            wv_f = np.matmul(AB,wv_i)
            t_c += dt
            wv_f[0] = ex_sln(x_i,l,t_c) # boundaries
            wv_f[-1] = ex_sln(x_f,l,t_c)
            wv_i = wv_f[:]

        # exact solution
        ex_wv = []
        for x_c in x_vals:
            v = ex_sln(x_c,l,t_f)
            ex_wv.append(v)

        # find norms
        l2norm = l2_norm(wv_f,ex_wv,h)
        mxnorm = max_norm(wv_f,ex_wv)

        # format for outputs
        output = "grph_h{}_m{}_lm{}.png".format(h,Mu,lmbda)
        title = "Crank-Nicolson \n h = {}, mu = {}, Lambda = {}".format(h,Mu,lmbda)
        print("h = {} | L2 Norm = {} | Supremum Norm = {} | mu = {} | lmba = {}".format(h,l2norm,mxnorm,Mu,lmbda))
        graph(output,title,x_vals,wv_f,ex_wv)

################################################################################
