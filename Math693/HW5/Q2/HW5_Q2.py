#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Program written for Homework 5, MATH 693B. Problem outline from
# Strikwerda 3.1.2

# Program solves the advection equation in 1D using the Forward-Time Backwards-
# Space (FTBS) scheme and the Lax-Wendroff (LW) scheme. The program will do this
# at a variety of values of h in order to compare order of accuracy against
# nominal values. Exact solution has been provided.

# Author: Jon Parsons
# Date: 4-15-20

################################################################################
#Functions
################################################################################
# Function containing the initial condition of ODE u, as well as exact solution
def u_ex(x):
    # x - point of interest

    sln = np.sin(np.pi*2*x)

    return sln

################################################################################
# Function containing the FTBS scheme
def FTBS(u_nm,u_nmp,k,h):
    # u_nm  - Value at local point
    # u_nmp - Value at point m-1
    # k     - Stepsize in time
    # h     - Stepsize in space

    u_n1 = u_nm - (k/h)*(u_nm - u_nmp)

    return u_n1

################################################################################
# Function containing the LW scheme
def LW(u_nm,u_nmp,u_nmx,k,h):
    # u_nm  - Value at local point
    # u_nmp - Value at point m-1
    # u_nmx - Value at point m+1
    # k     - Stepsize in time
    # h     - Stepsize in space

    u_n1 = u_nm - 0.5*(k/h)*(u_nmx-u_nmp) + 0.5*(k*k)/(h*h)*(u_nmx-2*u_nm+u_nmp)

    return u_n1

################################################################################
# Function contains code to find the L2 norm
def L2_norm(v,h):
    # v - Vector containing the numerical solution
    # e - Vector containing the exact solution
    # h - stepsize

    sum = 0
    for i in range(len(v)-1):
        sum += (v[i]*v[i]) # - e[i])**2

    l2 = np.sqrt(sum*h)

    return l2

################################################################################
# Function contains code to find the Max norm
def MX_norm(v):
    # v - Vector containing the numerical solution
    # e - Vector containing the exact solution

    diff = 0.0
    for i in range(len(v)-1):
        sub = abs(v[i])# - e[i]
        if sub > diff:
            diff = sub

    mx = abs(diff)

    return mx

################################################################################
# Function contains code to find the order of accuracy. Arguements must be of
# equal length
def acc_ord(v):
    # err - vector containing the error from the norms

    sum = 0
    N = 1
    while N < len(v)-1:
        ord_c = np.log2((v[N-1] - v[N])/(v[N] - v[N+1]))
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
    # y_vals  - holds FTBS calculated values
    # y2_vals - holds LW calculated values
    # y3_vals - holds exact values

    # X range limits
    x_l = min(x_vals)
    x_h = max(x_vals)

    # y range limits. Will auto-adjust based on lowest value of the exact sln
    y_l = min(y3_vals)
    y_h = max(y3_vals)


    plt.plot(x_vals,y_vals)
    plt.plot(x_vals,y2_vals)
    plt.plot(x_vals,y3_vals)
    plt.legend(["FTBS","LW","Exact"], loc='upper right')
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
h_vals = [1/10,1/20,1/40,1/80]
lmda = 0.8

# Holds values to print at if you dont want them all
t_print = [40,80,120]

# Initialize values
x_i = -1.0 # minimum x
x_f = 1.0 # maximum x
t_i = 0.0 # minimum time
t_f = 1.2 # maximum time, extended due to gif maker

lw_l2_norms = []
lw_mx_norms = []
ftbs_l2_norms = []
ftbs_mx_norms = []

q = 0
for h in h_vals:
    q += 1
    k = h*lmda
    print(h,k,lmda)

    x_vals = []
    ftbs_p = []
    lw_p = []
    ex_n = []
    ex_p = []
    ftbs_n = []
    lw_n = []

    # fill initial conditions
    x_c = x_i
    while x_c <= x_f:
        x_vals.append(x_c)
        ex_p.append(u_ex(x_c))
        ftbs_p.append(u_ex(x_c))
        lw_p.append(u_ex(x_c))
        x_c += h

    N = len(x_vals)

    t_c = t_i
    r = 0
    output = "grph_t"+str(r).zfill(3)+str(q)+".jpg".format(r,q)
    title = "Advecting Wave \n t = %.3f"% t_c
    graph(output,title,x_vals,ftbs_p,lw_p,ex_p)

    # Update system
    r += 1
    while t_c <= t_f + k:

        # Boundary conditions on the left
        lw_n.append(lw_p[-1])
        ftbs_n.append(ftbs_p[-1])

        # Update through domain
        for i in range(1,N-1):
            lw_n.append(LW(lw_p[i],lw_p[i-1],lw_p[i+1],k,h))
            ftbs_n.append(FTBS(ftbs_p[i],ftbs_p[i-1],k,h))
        for valx in x_vals:
            ex_n.append(u_ex(valx-t_c))

        # Boundary conditions on right
        lw_n.append(LW(lw_n[0],lw_n[-1],lw_n[1],k,h))
        ftbs_n.append(FTBS(ftbs_n[0],ftbs_n[-1],k,h))

        # format for output graphs
        output = "grph_t"+str(r).zfill(3)+"h"+str(q)+".jpg".format(r,q)
        title = "Advecting Wave \n t = %.3f"% t_c
        if  h == 1/80:
            graph(output,title,x_vals,ftbs_n,lw_n,ex_n)

        # Reset vectors
        lw_p = lw_n[:]
        ftbs_p = ftbs_n[:]
        ex_p = ex_n[:]

        lw_n.clear()
        ftbs_n.clear()
        ex_n.clear()

        t_c += k
        r += 1

        if r == round(1/h):
            lw_l2_norms.append(L2_norm(lw_p,h))
            lw_mx_norms.append(MX_norm(lw_p))
            ftbs_l2_norms.append(L2_norm(ftbs_p,h))
            ftbs_mx_norms.append(MX_norm(ftbs_p))


# Creates a gif of the outputs. comment out if unwanted.
#os.system("ffmpeg -y -i 'grph_t%03dh1.jpg' coupled_eqns1.gif")
#os.system("ffmpeg -y -i 'grph_t%03dh2.jpg' coupled_eqns2.gif")
#os.system("ffmpeg -y -i 'grph_t%03dh3.jpg' coupled_eqns3.gif")
#os.system("ffmpeg -y -i 'grph_t%03dh4.jpg' coupled_eqns4.gif")

#os.system("rm -f *.jpg")


print(ftbs_l2_norms)
print(ftbs_mx_norms)
print(lw_l2_norms)
print(lw_mx_norms)

print("Finding Error in FTBS L2")
ftbs_l2_ord = acc_ord(ftbs_l2_norms)
print("Finding Error in FTBS MAX")
ftbs_mx_ord = acc_ord(ftbs_mx_norms)
print("Finding Error in LW L2")
lw_l2_ord = acc_ord(lw_l2_norms)
print("Finding Error in LW MAX")
lw_mx_ord = acc_ord(lw_mx_norms)

print("FTBS Order of Accuracy from L2: ",ftbs_l2_ord)
print("FTBS Order of Accuracy from Max: ",ftbs_mx_ord)
print("LW Order of Accuracy from L2: ",lw_l2_ord)
print("LW Order of Accuracy from Max: ",lw_mx_ord)

################################################################################
