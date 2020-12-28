#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Program written for Final Project, Math 693b.
# Tests succesful implementation of the schemes to solve the hyperbolic
# wave equation

# This version implements a smoothly changing wave speed in space.

# Author: Jon Parsons
# Date: 4-22-20

################################################################################
#Functions
################################################################################
# Function containing the exact solution and the Boundary conditions. Solves
# for each point individually for the parabolic schemes
def u_ex(x,t):
    # x - point of interest
    # t - current time

    Pi = math.pi
    sum = 0
    for i in range(100):
        l_md = 2*i + 1
        sum += ((-1)**i)*(math.cos(Pi*l_md*x)/(Pi*(2*i+1)))*math.exp(-((Pi*l_md)**2)*t)

    sln = 0.5 + 2*sum
    return sln
################################################################################
################################################################################
# Function containing the exact solution and the Boundary conditions. Solves
# for each point individually for the parabolic schemes
def ini(x,t):
    # x - point of interest

    if abs(x) < 0.5:
        sln = 1
    elif abs(x) == 0.5:
        sln = 0.5
    elif abs(x) > 0.5:
        sln = 0

    return sln

################################################################################
# Function containing the FTFS scheme
def FTCS(u_nm,u_nmp,u_nmx,a,mu):
    # u_nm  - local point
    # u_nmp - left of point
    # u_nmx - right of point
    # mu    - value of mu, k/h^2
    # a     - local wavespeed

    u_n1 = u_nm + a*mu*(u_nmx - 2*u_nm + u_nmp)

    return u_n1

################################################################################
# Function containing the CN matrix fills - 1d version
def CN_mat(N,a,mu):
    # N   - Size of matrix
    # a   - wavespeed
    # mu  - k/h^2

    n = round(N)

    A = np.zeros((n,n))
    B = np.zeros((n,n))

    mu = 0.5*mu

    i = 0
    while i < n:
        if i == 0:
            A[i,i] = 1
            B[i,i] = 1
        elif i == n-1:
            A[i,i] = 1
            B[i,i] = 1
        else:
            A[i,i-1] = -mu
            A[i,i]   = 1 + 2*mu
            A[i,i+1] = -mu
            B[i,i-1] = mu
            B[i,i]   = 1 - 2*mu
            B[i,i+1] = mu
        i += 1

    return A, B

################################################################################
# Function contains code to find the L2 norm
def L2_norm(v,e,h):
    # v - Vector containing the numerical solution
    # e - Vector containing the exact solution
    # h - stepsize

    sum = 0
    for i in range(len(v)-1):
        sum += (e[i] - v[i])**2

    l2 = np.sqrt(sum*h)

    return l2

################################################################################
# Function contains code to find the order of accuracy. Arguements must be of
# equal length
def acc_ord(v):
    # v - vector containing the error from the norms

    sum = 0
    N = 1
    while N < len(v):
        ord_c = np.log(v[N-1]/v[N])/np.log(2)
        sum += ord_c
        print(ord_c)
        N += 1

    ord = sum/(N-1)

    return ord


################################################################################
# Graphing function, handles graphing of schemes with exact solution
def graph(out,title,x1_vals,y_vals,y2_vals):
    # out    - title of output file
    # title  - graph title
    # x_vals - holds x axis values
    # y_vals - holds numerical values
    # y2_vals- holds exact values

    # X range limits
    x_l = min(x_vals)
    x_h = max(x_vals)
    y_l = -1.5
    y_h = 1.5

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

print("Hello. Solving Systems")

# Common variables
h_vals = [1/10,1/20,1/40]

ftfs_l2_norms = np.zeros(len(h_vals))
cn_l2_norms = np.zeros(len(h_vals))

# Initialize values
x_i = -1.0 # minimum x
x_f = 1.0 # maximum x
t_i = 0.0 # minimum time
t_f = 0.1 # maximum time

q = 0
for h in h_vals:
    q += 1
    mu = 0.5
    k = h*h*mu
    print(h,k,mu)

    N = round((x_f-x_i)/h)

    x_vals = np.zeros(N+1)
    u_ftfs_p = np.zeros(N+1)
    u_ftfs_n = np.zeros(N+1)

    A, B = CN_mat(N+1,1,mu)

    u_cn_p = np.zeros(N+1)
    u_cn_i = np.zeros(N+1)
    u_cn_n = np.zeros(N+1)
    u_exact = np.zeros(N+1)

    # fill initial conditions
    x_c = x_i
    t_c = t_i

    i = 0
    while i <= N:
        x_vals[i] = x_i + i*h
        u_ftfs_p[i] = ini(x_vals[i],t_c)
        u_cn_p[i] = ini(x_vals[i],t_c)
        i += 1

    # Update system
    t_c = t_i + k
    r = 0
    while t_c <= t_f + k:

        u_ftfs_n[0] = u_ex(x_i,t_c)
        u_exact[0] = u_ex(x_i,t_c)

        i = 1
        while i < N:
            u_ftfs_n[i] = FTCS(u_ftfs_p[i],u_ftfs_p[i-1],u_ftfs_p[i+1],1,mu)
            u_exact[i] = u_ex(x_vals[i],t_c)
            i += 1

        u_ftfs_n[N] = u_ex(x_f,t_c)
        u_exact[-1] = u_ex(x_f,t_c)

        output = "ftcs_t"+str(r).zfill(3)+"_h"+str(q).zfill(1)+ \
                                                              ".jpg".format(r,q)
        title = "FTCS \n t = %.3f"% t_c

        graph(output,title,x_vals,u_ftfs_n,u_exact)
        # CN Part
        u_cn_i = np.matmul(B,u_cn_p)

        u_cn_n = np.linalg.solve(A,u_cn_i)

        u_cn_n[0] = u_ex(x_i,t_c)
        u_cn_n[-1] = u_ex(x_f,t_c)

        output = "cn_t"+str(r).zfill(3)+"_h"+str(q).zfill(1)+".jpg".format(r,q)
        title = "CN \n t = %.3f"% t_c

        graph(output,title,x_vals,u_cn_n,u_exact)

        for l, val in enumerate(u_ftfs_n):
            u_ftfs_p[l] = val
            u_cn_p[l] = u_cn_n[l]

        u_cn_p[0] = u_ex(x_i,t_c)
        u_cn_p[-1] = u_ex(x_f,t_c)

        t_c += k
        r += 1

    ftfs_l2_norms[q-1] = L2_norm(u_ftfs_p,u_exact,h)
    cn_l2_norms[q-1]   = L2_norm(u_cn_p,u_exact,h)

ftfs_ord = acc_ord(ftfs_l2_norms)
cn_ord = acc_ord(cn_l2_norms)

print("FTCS",ftfs_ord)
print("FTCS L2 Norms: ",ftfs_l2_norms)
print("CN",cn_ord)
print("CN L2 Norms: ", cn_l2_norms)
# Creates a gif of the outputs. comment out if unwanted.
os.system("ffmpeg -y -i 'ftcs_t%03d_h1.jpg' ftcs_h1.gif")
os.system("ffmpeg -y -i 'ftcs_t%03d_h2.jpg' ftcs_h2.gif")
os.system("ffmpeg -y -i 'ftcs_t%03d_h3.jpg' ftcs_h3.gif")

os.system("ffmpeg -y -i 'cn_t%03d_h1.jpg' cn_h1.gif")
os.system("ffmpeg -y -i 'cn_t%03d_h2.jpg' cn_h2.gif")
os.system("ffmpeg -y -i 'cn_t%03d_h3.jpg' cn_h3.gif")

os.system("rm -f *.jpg")
################################################################################
