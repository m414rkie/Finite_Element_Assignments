#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Program written for Homework 1, MATH 693B. Problem outline from
# Strikwerda 1.3.1

# Program will calculate the evolution of a simple wave using four schemes:
# 1. Forward-Time Backwards-Space (FTBS)
# 2. Forward-Time Central-Space (FTCS)
# 3. Lax-Friedrichs (LFx)
# 4. Leapfrog (LFg)

# The schemes will be run with parameter lambda at lm = 0.8 or 1.6. Outputs will
# be graphed at the final timestep.

# Author: Jon Parsons
# Date: 2-11-20

################################################################################
#Functions
################################################################################
# Function containing the initial conditions
def f(x):
    # x: spatial coordinate
    pi = np.pi

    f_x = np.cos(pi*x)**2
    return f_x

################################################################################
# Function containing the FTBS algorithm, wave velocity is 1
def FTBS(u_m,u_mp,lm):
    # u_m: current time spatial coordinate of interest
    # u_mp: current time spatial coordinate immediately right of u_m
    # lm: lambda, dt/dx

    u_f = u_m - lm*(u_m - u_mp)
    return u_f

################################################################################
# Function containing the FTCS algorithm, wave velocity is 1
def FTCS(u_m,u_mn,u_mp,lm):
    # u_m: current time spatial coordinate of interest
    # u_mn: current time spatial coordinate immediately left of u_m
    # u_mp: current time spatial coordinate immediately right of u_m
    # lm: lambda, dt/dx

    u_f = u_m - 0.5*lm*(u_mn - u_mp)
    return u_f

################################################################################
# Function containing the LFx algorithm, wave velocity is 1
def LFx(u_m,u_mn,u_mp,lm):
    # u_m: current time spatial coordinate of interest
    # u_mn: current time spatial coordinate immediately left of u_m
    # u_mp: current time spatial coordinate immediately right of u_m
    # lm: lambda, dt/dx

    u_f = 0.5*(u_mn + u_mp) - 0.5*lm*(u_mn - u_mp)
    return u_f

################################################################################
# Function containing the LFg algorithm, wave velocity is 1
def LFg(u_np,u_mn,u_mp,lm):
    # u_np: previous time spatial coordinate of interest
    # u_mn: current time spatial coordinate immediately left of u_m
    # u_mp: current time spatial coordinate immediately right of u_m
    # lm: lambda, dt/dx

    u_f = u_np - lm*(u_mn - u_mp)
    return u_f

################################################################################
# Graphing function, handles graphing of schemes
def graph(method,lm,h,x_vals,y_vals):

    title = "{} Method \n H = 1/{}, lambda = {}".format(method,h,lm)
    name = "{}_{}_{}.png".format(method,h,lm)

    plt.plot(x_vals,y_vals)
    plt.xlim((-1.0,3.0))
    plt.ylim((0.0,1.0))
    plt.ylabel('Wave Magnitude')
    plt.xlabel('X')
    plt.title(title)
    plt.savefig(name)
    plt.clf()
################################################################################
# Main - shifts between methods, changes h changes lambda when appropriate

print("Hello. Making Waves.")
# Common variables
h = [10,20,40]
x_i = -1.0
x_f = 3.0
t_i = 0.0
t_f = 2.4

x_vals = []

# FTBS Section
print("FTBS")

for item in h:
    # Initialize Variables
    lmda = 0.8
    grdx = (x_f - x_i)/(1/item)
    grdt = lmda*grdx

    wv_i = []

    print("Gridsize: {}, {}".format(grdx,grdt))
    # Fill initial values
    for i in range(int(grdx+1)):
        x_c = x_i + i*(1/item)
        x_vals.append(x_c)
        if abs(x_c) <= 0.5:
            wv = f(x_c)
        else:
            wv = 0.0

        wv_i.append(wv)

    wv_f = [0.0]
    print("Iterating Time")
    for t in range(int(grdt+1)):
        for j in range(int(grdx+1)):
            wv_n = FTBS(wv_i[int(j)],wv_i[int(j)-1],lmda)
            wv_f.append(wv_n)

        wv_i = [0.0]
        wv_i = wv_f[:]
        wv_i[0] = 0.0 # Enforce Boundary Conditions
        wv_f = []

    graph('FTBS',lmda,item,x_vals,wv_i)

    wv_f.clear()
    wv_i.clear()
    x_vals.clear()

x_vals = []

################################################################################
# FTCS Section
print("FTCS")

for item in h:
    # Initialize variables
    lmda = 0.1
    grdx = (x_f - x_i)/(1/item)
    grdt = lmda*grdx

    wv_i = []

    print("Gridsize: {}, {}".format(grdx,grdt))
    # Fill initial conditions
    for i in range(int(grdx+1)):
        x_c = x_i + i*(1/item)
        x_vals.append(x_c)
        if abs(x_c) <= 0.5:
            wv = f(x_c)
        else:
            wv = 0.0

        wv_i.append(wv)

    wv_f = [0.0]
    print("Iterating Time")
    for t in range(int(grdt+1)):
        for j in range(int(grdx+1)):
            if j < int(grdx):
                wv_n = FTCS(wv_i[j],wv_i[j+1],wv_i[j-1],lmda)
            else:
                wv_n = wv_f[-1]

            wv_f.append(wv_n)

        wv_i = [0.0]
        wv_i = wv_f[:]
        wv_i[0] = 0.0 # Enforce Boundary Conditions
        wv_f = []

    graph('FTCS',lmda,item,x_vals,wv_i)

    wv_f.clear()
    wv_i.clear()
    x_vals.clear()


################################################################################
# LFx Section
print("LFx")
# This method gets two values for lambda
lm_lst = [0.8,1.6]
for item in lm_lst:
    lmda = item

    for item in h:
        grdx = (x_f - x_i)/(1/item)
        grdt = lmda*grdx

        wv_i = []

        print("Gridsize: {}, {}".format(grdx,grdt))
        # Fill initial conditions
        for i in range(int(grdx+1)):
            x_c = x_i + i*(1/item)
            x_vals.append(x_c)
            if abs(x_c) <= 0.5:
                wv = f(x_c)
            else:
                wv = 0.0

            wv_i.append(wv)

        wv_f = [0.0]
        print("Iterating Time")
        for t in range(int(grdt+1)):
            for j in range(int(grdx+1)):
                if j < int(grdx):
                    wv_n = LFx(wv_i[j],wv_i[int(j+1)],wv_i[int(j)-1],lmda)
                else:
                    wv_n = wv_f[-1]#LFx(wv_i[j],wv_i[int(j)],wv_i[int(j)-1],lmda)

                wv_f.append(wv_n)

            wv_i = [0.0]
            wv_i = wv_f[:]
            wv_i[0] = 0.0 # Enforce Boundary Conditions
            wv_f = []

        graph('Lfx',lmda,item,x_vals,wv_i)

        wv_f.clear()
        wv_i.clear()
        x_vals.clear()


################################################################################
# LFg Section
print("LFg")

for item in h:
    # Initialize values
    lmda = 0.8
    grdx = (x_f - x_i)/(1/item)
    grdt = lmda*grdx

    print("Gridsize: {}, {}".format(grdx,grdt))

    wv_p = []
    # Fill initial conditions
    for i in range(int(grdx+1)):
        x_c = x_i + i*(1/item)
        x_vals.append(x_c)
        if abs(x_c) <= 0.5:
            wv = f(x_c)
        else:
            wv = 0.0

        wv_p.append(wv)

    # Iterate first timestep using FTCS
    wv_f = [0.0]
    wv_i = [0.0]
    for j in range(int(grdx+1)):
        if j < int(grdx):
            wv_n = FTCS(wv_p[j],wv_p[j+1],wv_p[j-1],lmda)
        else:
            wv_n = wv_p[-1]

        wv_i.append(wv_n)

    print("Iterating Time")
    for t in range(int(grdt)):
        for j in range(int(grdx+1)):
            if j < int(grdx):
                wv_n = LFg(wv_p[j],wv_i[int(j+1)],wv_i[int(j)-1],lmda)
            else:
                wv_n = LFg(wv_p[j],wv_i[int(j)],wv_i[int(j)-1],lmda)

            wv_f.append(wv_n)

        wv_p = wv_i[:]
        wv_p[0] = 0.0 # Enforce Boundary Conditions
        wv_i = wv_f[:]
        wv_i[0] = 0.0 # Enforce Boundary Conditions
        wv_f = []

    graph('Lfg',lmda,item,x_vals,wv_i)

    wv_f.clear()
    wv_i.clear()
    wv_p.clear()
    x_vals.clear()
