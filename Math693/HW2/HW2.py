#!usr/bin/python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Program written for Homework 2, MATH 693B. Problem outline from
# Strikwerda 3.4.1

# Program will provide solutions for the 1-way wave equation using 4 boundary
# conditions.
# A) V^{n+1}_{M} = 2V^{n+1}_{M-1} - V^{n+1}_{M-2}
# B) V^{n+1}_{M} = 0
# C) V^{n+1}_{M} = V^{n}_{M-1} and V^{n+1}_{0} = 2V^{n+1}_{1} - V^{n+1}_{2}
# D) V^{n+1}_{M} = 2V^{n}_{M-1} - V^{n-1}_{M-2}

# The schemes will be run with parameter lambda at lm = 0.5. Outputs will
# be graphed at the final timestep. The initial values will be u(0,0) = 1.
# Region of interest: x in [0,1], t in [0,1.0]

# Author: Jon Parsons
# Date: 2-25-20

################################################################################
#Functions
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
    plt.xlim((0.0,1.0))
    plt.ylim((0.0,2.0))
    plt.ylabel('Wave Magnitude')
    plt.xlabel('X')
    plt.title(title)
    plt.savefig(name)
    plt.clf()

################################################################################
# Main
################################################################################

print("Hello. Making Waves.")
# Common variables
h = 50
x_i = 0.0
x_f = 1.0
t_i = 0.0
t_f = 1.0

x_vals = []

# Initialize values
lmda = 0.3
grdx = (x_f - x_i)/(1/h)
grdt = (t_f - t_i)/(lmda/(h))

print("Gridsize: {}, {}".format(grdx,grdt))

# Part A - Boundary Condition 1
wv_p = []
# Fill initial conditions
for i in range(int(grdx)):
    x_c = x_i + i*(1/h)
    x_vals.append(x_c)
    if i == 1:
        wv = 1.0
    else:
        wv = 0.0

    wv_p.append(wv)

# First iteration will take values from u(0,x) = u(t_1,x)
wv_f = []
wv_i = wv_p[:]

print("Iterating Time Part A")
for t in range(int(grdt)):
    for j in range(int(grdx)):
        if j < int(grdx-1):
            wv_n = LFg(wv_p[j],wv_i[int(j+1)],wv_i[int(j-1)],lmda)
        else:
            wv_n = 2.0*wv_f[-1] - wv_f[-2]

        wv_f.append(wv_n)

    wv_p = wv_i[:]
    wv_p[0] = 1.0 # Enforce Boundary Conditions
    wv_i = wv_f[:]
    wv_i[0] = 1.0 # Enforce Boundary Conditions
################################################################################

    wv_f = []

graph('Part_A',lmda,h,x_vals,wv_i)

wv_f.clear()
wv_i.clear()
wv_p.clear()

################################################################################
################################################################################

# Part B - Boundary Condition 1
wv_p = []
# Fill initial conditions
for i in range(int(grdx)):
    if i == 1:
        wv = 1.0
    else:
        wv = 0.0

    wv_p.append(wv)

# First iteration will take values from u(0,x) = u(t_1,x)
wv_f = []
wv_i = wv_p[:]

print("Iterating Time Part B")
for t in range(int(grdt)):
    for j in range(int(grdx)):
        if j < int(grdx-1):
            wv_n = LFg(wv_p[j],wv_i[int(j+1)],wv_i[int(j-1)],lmda)
        else:
            wv_n = 0.0

        wv_f.append(wv_n)

    wv_p = wv_i[:]
    wv_p[0] = 1.0 # Enforce Boundary Conditions
    wv_p[-1] = 0.0 # Enforce Boundary Conditions
    wv_i = wv_f[:]
    wv_i[0] = 1.0 # Enforce Boundary Conditions
    wv_i[-1] = 0.0 # Enforce Boundary Conditions

################################################################################

    wv_f = []

graph('Part_B',lmda,h,x_vals,wv_i)

wv_f.clear()
wv_i.clear()
wv_p.clear()

################################################################################
################################################################################

# Part C - Boundary Condition 1
wv_p = []
# Fill initial conditions
for i in range(int(grdx)):
    if i == 1:
        wv = 1.0
    else:
        wv = 0.0

    wv_p.append(wv)

# First iteration will take values from u(0,x) = u(t_1,x)
wv_f = []
wv_i = wv_p[:]

print("Iterating Time Part C")
for t in range(int(grdt)):
    for j in range(int(grdx-1)):
        if j < int(grdx-2):
            wv_n = LFg(wv_p[j],wv_i[int(j+1)],wv_i[int(j-1)],lmda)
        elif j == grdx-1:
            wv_n = wv_i[-2]

        wv_f.append(wv_n)

    bnd_o = 2.0*wv_f[0] - wv_f[1]
    wv_f.insert(0,bnd_o)

    wv_p = wv_i[:]
    wv_i = wv_f[:]

################################################################################

    wv_f = []

graph('Part_C',lmda,h,x_vals,wv_i)

wv_f.clear()
wv_i.clear()
wv_p.clear()

################################################################################
################################################################################

# Part D - Boundary Condition 1
wv_p = []
# Fill initial conditions
for i in range(int(grdx)):
    if i == 1:
        wv = 1.0
    else:
        wv = 0.0

    wv_p.append(wv)

# First iteration will take values from u(0,x) = u(t_1,x)
wv_f = []
wv_i = wv_p[:]
bnd_mp = 0.0

print("Iterating Time Part D")
for t in range(int(grdt)):
    for j in range(int(grdx)):
        if j < int(grdx-1):
            wv_n = LFg(wv_p[j],wv_i[int(j+1)],wv_i[int(j-1)],lmda)
        else:
            wv_n = bnd_mp

        wv_f.append(wv_n)
    bnd_mp = wv_f[-2]

    wv_p = wv_i[:]
    wv_p[0] = 1.0
    wv_i = wv_f[:]
    wv_i[0] = 1.0
################################################################################

    wv_f = []

graph('Part_D',lmda,h,x_vals,wv_i)

wv_f.clear()
wv_i.clear()
wv_p.clear()
x_vals.clear()
