#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Program written for Homework 6, MATH 693B. Problem outline from
# Strikwerda 3.5.2

# Program solves the equation u_t + u_x + u = 0 using the Crank-Nicolson
# method. 2 different initial conditions and 2 different boundary conditions
# will be used and compared. The system will be solved using the Thomas
# algorithm

# Author: Jon Parsons
# Date: 4-20-20

################################################################################
#Functions
################################################################################
# Function containing the initial condition of 1st IC, as well as exact solution
def ex_1(x,t):
    # x - point of interest

    sln = np.exp(-t)*np.sin(np.pi*(x-t))

    return sln

################################################################################
# Function containing the initial condition of 2nd IC, as well as exact solution
def ex_2(x,t):
    # x - point of interest

    sln = np.exp(-t)*np.cos(np.pi*(x-t))

    if 0 > sln:
        return 0
    else:
        return sln

################################################################################
# Function that fills the A and B arrays of the Crank-Nicolson method, versions
# for both boundary conditions will be returned. LHS matrices are reduced due to
# inversion method

def mat_fill(k,h,n):
    # k  - Stepsize in time
    # h  - Stepsize in space
    # n  - number of elements

    A_a = np.zeros(round(n)-1)
    A_b = np.zeros(round(n))
    A_c = np.zeros(round(n)-1)

    B = np.zeros((round(n),round(n)))

    for i in range(round(n)):
        if i == 0:
            A_b[i] = 1
            B[i,i] = 1
        elif i == n-1:
            A_b[i] = 1
            B[i,i] = 1
        else:
            A_a[i] = -k/(4*h)
            A_b[i] = 1
            A_c[i] = k/(4*h)

            B[i,i-1] = k/(4*h)
            B[i,i]   = 1 - k
            B[i,i+1] = -k/(4*h)

    A_a[0] = -k/(4*h)
    A_a[-1] = 0
    A_c[-1] = k/(4*h)
    A_c[0] = 0

    return A_a, A_b, A_c, B

################################################################################
# Function to initialize the waves
def ini(dx,x_i,t_i,n):
    # dx  - Stepsize in x
    # x_i - Initial value of x
    # t_i - Initial value in time
    # n   - Number of gridpoints

    u_ic1_b1 = np.zeros(n)
    u_ic1_b2 = np.zeros(n)
    u_ic2_b1 = np.zeros(n)
    u_ic2_b2 = np.zeros(n)
    ex_ic1   = np.zeros(n)
    ex_ic2   = np.zeros(n)
    x_vals   = np.zeros(n)

    x_c = x_i
    for i in range(n):
        x_vals[i] = x_c
        ic1 = ex_1(x_c,t_i)
        ic2 = ex_2(x_c,t_i)
        u_ic1_b1[i] = ic1
        u_ic1_b2[i] = ic1
        u_ic2_b1[i] = ic2
        u_ic2_b2[i] = ic2
        ex_ic1[i] = ic1
        ex_ic2[i] = ic2
        x_c += dx


    return u_ic1_b1, u_ic1_b2, u_ic2_b1, u_ic2_b2, ex_ic1, ex_ic2, x_vals

################################################################################
# Function containing the Thompson algrothm
def Thomas(a,b,c,d):
    # a - Vector containing the sub-diagonal values
    # b - Vector containing the main diagonal values
    # c - Vector containing the super-diagonal values
    # d - Vector containing the RHS values

    n = len(d)
    cp = np.zeros(n-1) # c prime
    dp = np.zeros(n) # d prime
    u_nv = np.zeros(n) # inversed vector

    cp[0] = (c[0]/b[0])
    dp[0] = (d[0]/b[0])

    for i in range(1,n-2):
        cp[i] = c[i]/(b[i] - a[i]*cp[i-1])
    for i in range(1,n-1):
        dp[i] = (d[i] - a[i]*dp[i-1])/(b[i] - a[i]*cp[i-1])

    u_nv[n-1] = dp[n-1]
    for i in range(n-2,-1,-1):
        u_nv[i] = dp[i] - cp[i]*u_nv[i+1]

    return u_nv

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
def graph(out,title,x1_vals,y_vals,y2_vals,y3_vals):
    # out     - title of output file
    # title   - graph title
    # x_vals  - holds x axis values
    # y_vals  - holds  BC 1
    # y2_vals - holds  BC 2
    # y3_vals - holds Exact

    # X range limits
    x_l = min(x_vals)
    x_h = max(x_vals)

    # y range limits. Will auto-adjust based on lowest value of the exact sln
    y_l = -1.1
    y_h = 1.1

    plt.plot(x_vals,y_vals)
    plt.plot(x_vals,y2_vals)
    plt.plot(x_vals,y3_vals)
    plt.legend(["BC 1","BC 2","Exact"], loc='upper right')
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

print("Hello. Moving waves.")

# Common variables
h_vals = [1/10,1/20,1/40]
lmda = 1.0

# Initialize values
x_i = -1.0 # minimum x
x_f = 1.0 # maximum x
t_i = 0.0 # minimum time
t_f = 1.0 # maximum time, extended due to gif maker

l2_norms_ic1_bc1 = np.zeros(len(h_vals))
l2_norms_ic1_bc2 = np.zeros(len(h_vals))
l2_norms_ic2_bc1 = np.zeros(len(h_vals))
l2_norms_ic2_bc2 = np.zeros(len(h_vals))

q = 0
for h in h_vals:
    q += 1
    k = h*lmda
    print("Solving for h = ",h)
    N = round((x_f-x_i)/h) + 1

    # Matrices
    A_a, A_b, A_c, B = mat_fill(k,h,N)
    # Initial Vectors
    u_ic1_b1, u_ic1_b2, u_ic2_b1, u_ic2_b2, ex_ic1, ex_ic2, x_vals = \
                                                                ini(h,x_i,t_i,N)

    t_c = t_i
    r = 0
    # Initial values
    output_ic1 = "grph1_t"+str(r).zfill(3)+str(q)+".jpg".format(r,q)
    output_ic2 = "grph2_t"+str(r).zfill(3)+str(q)+".jpg".format(r,q)
    title_ic1 = "Initial Condition Set 1 \n t = %.3f"% t_c
    title_ic2 = "Initial Condition Set 2 \n t = %.3f"% t_c
    graph(output_ic1,title_ic1,x_vals,u_ic1_b1,u_ic1_b2,ex_ic1)
    graph(output_ic2,title_ic2,x_vals,u_ic2_b1,u_ic2_b2,ex_ic2)

    # Update system
    r += 1
    t_c += k
    while t_c <= t_f + k:
        # Initial Condition 1
        # BC 1
        u_ic1_b1_i = np.matmul(B,u_ic1_b1)
        u_ic1_b1_n = Thomas(A_a,A_b,A_c,u_ic1_b1_i)
        # BC 2
        u_ic1_b2_i = np.matmul(B,u_ic1_b2)
        u_ic1_b2_n = Thomas(A_a,A_b,A_c,u_ic1_b2_i)
        # Initial Condition 2
        # BC 1
        u_ic2_b1_i = np.matmul(B,u_ic2_b1)
        u_ic2_b1_n = Thomas(A_a,A_b,A_c,u_ic2_b1_i)
        # BC 2
        u_ic2_b2_i = np.matmul(B,u_ic2_b2)
        u_ic2_b2_n = Thomas(A_a,A_b,A_c,u_ic2_b2_i)

        # BC's for IC 1
        u_ic1_b1_n[N-1] = (u_ic1_b1[N-1]+lmda*u_ic1_b1_n[N-2])/(1+lmda+k)
        u_ic1_b2_n[N-1] = 2*u_ic1_b2_n[N-2] - u_ic1_b2_n[N-3]
        u_ic1_b1_n[0] = ex_1(x_i,t_c)
        u_ic1_b2_n[0] = ex_1(x_i,t_c)

        # BC's for IC 2
        u_ic2_b1_n[N-1] = (u_ic2_b1[N-1]+lmda*u_ic2_b1_n[N-2])/(1+lmda+k)
        u_ic2_b2_n[N-1] = 2*u_ic2_b2_n[N-2] - u_ic2_b2_n[N-3]
        u_ic2_b1_n[0] = ex_2(x_i,t_c)
        u_ic2_b2_n[0] = ex_2(x_i,t_c)

        for j, xx in enumerate(x_vals):
            ex_ic1[j] = ex_1(xx,t_c)
            ex_ic2[j] = ex_2(xx,t_c)

        output_ic1 = "grph1_t"+str(r).zfill(3)+"h"+str(q)+".jpg".format(r,q)
        output_ic2 = "grph2_t"+str(r).zfill(3)+"h"+str(q)+".jpg".format(r,q)
        title_ic1 = "Initial Condition Set 1 \n t = %.3f"% t_c
        title_ic2 = "Initial Condition Set 2 \n t = %.3f"% t_c
        graph(output_ic1,title_ic1,x_vals,u_ic1_b1_n,u_ic1_b2_n,ex_ic1)
        graph(output_ic2,title_ic2,x_vals,u_ic2_b1_n,u_ic2_b2_n,ex_ic2)

        # Update vectors
        for i, val in enumerate(u_ic1_b1_n):
            u_ic1_b1[i] = val
            u_ic1_b2[i] = u_ic1_b2_n[i]
            u_ic2_b1[i] = u_ic2_b1_n[i]
            u_ic2_b2[i] = u_ic2_b2_n[i]

        t_c += k
        r += 1

    l2_norms_ic1_bc1[q-1] = (L2_norm(u_ic1_b1,ex_ic1,h))
    l2_norms_ic1_bc2[q-1] = (L2_norm(u_ic1_b2,ex_ic1,h))
    l2_norms_ic2_bc1[q-1] = (L2_norm(u_ic2_b1,ex_ic2,h))
    l2_norms_ic2_bc2[q-1] = (L2_norm(u_ic2_b2,ex_ic2,h))

print("Determining Order of Accuracy")

ic1_B1_ord = acc_ord(l2_norms_ic1_bc1)
ic1_B2_ord = acc_ord(l2_norms_ic1_bc2)
ic2_B1_ord = acc_ord(l2_norms_ic2_bc1)
ic2_B2_ord = acc_ord(l2_norms_ic2_bc2)

print(l2_norms_ic1_bc1)
print(l2_norms_ic1_bc2)
print(l2_norms_ic2_bc1)
print(l2_norms_ic2_bc2)

print("Order for IC1, B1: ", ic1_B1_ord)
print("Order for IC1, B2: ", ic1_B2_ord)
print("Order for IC2, B1: ", ic2_B1_ord)
print("Order for IC2, B2: ", ic2_B2_ord)

# Creates a gif of the outputs. comment out if unwanted.
os.system("ffmpeg -y -t 3 -i 'grph1_t%03dh1.jpg' cn_eqns11.gif")
os.system("ffmpeg -y -t 3 -i 'grph2_t%03dh1.jpg' cn_eqns21.gif")
os.system("ffmpeg -y -t 3 -i 'grph1_t%03dh2.jpg' cn_eqns12.gif")
os.system("ffmpeg -y -t 3 -i 'grph2_t%03dh2.jpg' cn_eqns22.gif")
os.system("ffmpeg -y -t 3 -i 'grph1_t%03dh3.jpg' cn_eqns13.gif")
os.system("ffmpeg -y -t 3 -i 'grph2_t%03dh3.jpg' cn_eqns23.gif")

#os.system("rm -f *.jpg")

################################################################################
