#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Program written for Final Project, Math 693b.
# Tests succesful implementation of the schemes to solve the second order wave
# equation

# Author: Jon Parsons
# Date: 4-22-20

################################################################################
#Functions
################################################################################
# Function containing the initial condition of ODE u
def u_ex(x,t):
    # x - point of interest
    # t - time of interest
    # s - wavespeed

    sln = np.cos(x+t) + np.cos(x-t)

    return sln

################################################################################
# Function containing the formula for wavespeed
def speed(x_i,x_f,n):
    # x_i - initial x position
    # x_f - final x position
    # n   - size of grid

    s = np.zeros(n)
    ds = 2.1/n

    i = 0
    while i < n:
        s[i] = (0.5 + i*ds)
        i += 1

    return s

################################################################################
# Function containing the FTFS scheme
def CTCS(u_nm,u_nmp,u_nmx,u_npm,a,lmda):
    # u_nm  - local point
    # u_nmp - left of point
    # u_nmx - right of point
    # u_npm - point at previous time
    # lmda  - value of lambda
    # a     - local wavespeed

    cc = ((a*a))*lmda*lmda
    u_n1 = 2*u_nm - u_npm + cc*(u_nmx - 2*u_nm + u_nmp)

    return u_n1

################################################################################
# Function containing the CN matrix fills - 1d version
def twostp_mat(N,a,lmda):
    # N    - Size of matrix
    # a    - wavespeed
    # lmda - k/h

    n = round(N)

    A = np.zeros((n,n))
    B = np.zeros((n,n))

    i = 0
    while i < n:
        if i == 0:
            A[i,i] = 1
            B[i,i] = 1
        elif i == n-1:
            A[i,i] = 1
            B[i,i] = 1
        else:
            cp = 0.5*((a[i+1]*a[i+1]))*lmda*lmda
            cm = 0.5*((a[i-1]*a[i-1]))*lmda*lmda
            A[i,i-1] = -cp
            A[i,i]   = 1 + 2*cp
            A[i,i+1] = -cp
            B[i,i-1] = cm
            B[i,i]   = -(1 + 2*cm)
            B[i,i+1] = cm
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
    y_l = -2.2
    y_h = 2.2

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

ctcs_l2_norms = np.zeros(len(h_vals))
two_l2_norms = np.zeros(len(h_vals))

# Initialize values
x_i = -1.0 # minimum x
x_f = 1.0 # maximum x
t_i = 0.0 # minimum time
t_f = 4.5 # maximum time

q = 0
for h in h_vals:
    q += 1

    lmda = 0.5
    k = h*lmda

    print(h,k,lmda)

    N = round((x_f-x_i)/h)
    a = speed(x_i,x_f,N+1)

    x_vals = np.zeros(N+1)
    u_ctcs_p = np.zeros(N+1)
    u_ctcs_c = np.zeros(N+1)
    u_ctcs_n = np.zeros(N+1)

    u_2s_p = np.zeros(N+1)
    u_2s_c = np.zeros(N+1)
    u_2s_n = np.zeros(N+1)
    u_2s_i = np.zeros(N+1)

    u_exact = np.zeros(N+1)

    A, B = twostp_mat(N+1,a,lmda)
    A_inv = np.linalg.inv(A)


    x_c = x_i
    t_c = t_i + k

    # fill initial conditions
    i = 0
    while i <= N:
        x_vals[i] = x_i + i*h
        u_ctcs_p[i] = u_ex(x_vals[i],t_i)
        u_ctcs_c[i] = u_ex(x_vals[i],t_c)
        u_2s_p[i] = u_ex(x_vals[i],t_i)
        u_2s_c[i] = u_ex(x_vals[i],t_c)
        i += 1

    # Update system
    t_c = t_i + k
    r = 0
    while t_c <= t_f + k:

        i = 1
        u_ctcs_n[0] = u_ex(x_i,t_c)
        u_exact[0] = u_ex(x_i,t_c)


        while i < N:
            ab = a[i]
            u_ctcs_n[i] = CTCS(u_ctcs_c[i],u_ctcs_c[i-1],u_ctcs_c[i+1],\
                                                          u_ctcs_p[i],ab,lmda)
            u_exact[i] = u_ex(x_vals[i],t_c)
            i += 1

        u_ctcs_n[N] = u_ex(x_f,t_c)
        u_exact[N] = u_ex(x_f,t_c)

        output = "ctcs_t"+str(r).zfill(3)+"_h"+str(q).zfill(1) \
                                                             +".jpg".format(r,q)
        title = "CTCS \n t = %.3f"% t_c

        graph(output,title,x_vals,u_ctcs_n,u_exact)
        # 2 step Part
        u_2s_i = np.matmul(B,u_2s_p)
        u_2s_i = u_2s_i + 2*u_2s_c
        u_2s_i[0] = u_ex(x_i,t_c)
        u_2s_i[-1] = u_ex(x_f,t_c)
        u_2s_n = np.matmul(A_inv,u_2s_i)

        u_2s_n[0] = u_ex(x_i,t_c)
        u_2s_n[-1] = u_ex(x_f,t_c)

        output = "2s_t"+str(r).zfill(3)+"_h"+str(q).zfill(1)+".jpg".format(r,q)
        title = "2S \n t = %.3f"% t_c

        graph(output,title,x_vals,u_2s_n,u_exact)

        for l, val in enumerate(u_ctcs_n):
            u_ctcs_p[l] = u_ctcs_c[l]
            u_ctcs_c[l] = val
            u_2s_p[l] = u_2s_c[l]
            u_2s_c[l] = u_2s_n[l]

        t_c += k
        r += 1

    ctcs_l2_norms[q-1] = L2_norm(u_ctcs_p,u_exact,h)
    two_l2_norms[q-1] = L2_norm(u_2s_p,u_exact,h)

print(ctcs_l2_norms)
print(two_l2_norms)

ctcs_ord = acc_ord(ctcs_l2_norms)
two_ord = acc_ord(two_l2_norms)

print(ctcs_ord)
print(two_ord)
# Creates a gif of the outputs. comment out if unwanted.
os.system("ffmpeg -y -i 'ctcs_t%03d_h1.jpg' ctcs_h1.gif")
os.system("ffmpeg -y -i 'ctcs_t%03d_h2.jpg' ctcs_h2.gif")
os.system("ffmpeg -y -i 'ctcs_t%03d_h3.jpg' ctcs_h3.gif")

os.system("ffmpeg -y -i '2s_t%03d_h1.jpg' 2s_h1.gif")
os.system("ffmpeg -y -i '2s_t%03d_h2.jpg' 2s_h2.gif")
os.system("ffmpeg -y -i '2s_t%03d_h3.jpg' 2s_h3.gif")

#os.system("rm -f *.jpg")
os.system("rm -f *h1.jpg")
os.system("rm -f *h2.jpg")

################################################################################
