#!/usr/bin/python3

from scipy.interpolate import griddata
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import math

# This code takes known output files and plots the data within.

# 1.1 - Replaces the repetitive graphing code with a function
# 1.2 - Polynomial fit of error implemented

# Python code for HW4, COMP670

# Jon Parsons
# 11-1-19

def gifb(sz):
# Function which takes the 'X' delimited values of the wave, plots and saves
# to a gif. This is not an ideal implementation, more work is required for
# speed up.

    # Initialize
    frames = []
    cur_frame = []
    cur_frame.append([])
    cur_frame.append([])
    cur_frame.append([])

    count = []
    x = []
    y = []
    va = []
    ims = []

    fig = plt.figure()

    file = "time_flow_b.dat"
    x_lab = "X (m)"
    y_lab = "Y (m)"
    title = "Advection-Diffusion \n Concentration of Bacteria (1/mL)"

    fr = open(file,'r')

    # Read in data
    for line in fr:
        lin = line.split()
        # Save frames
        if lin[0] == "X":
            cur_frame[0] = x
            cur_frame[1] = y
            cur_frame[2] = va
            frames.append(cur_frame)
            cur_frame = []
            cur_frame.append([])
            cur_frame.append([])
            cur_frame.append([])
            x = []
            y = []
            va = []
        else: # Build frames
            x.append(float(lin[0]))
            y.append(float(lin[1]))
            va.append(float(lin[2]))


    # Create the plots
    for k,frm in enumerate(frames):

        plt.title(title)
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)
        im = plt.tricontourf(frm[0],frm[1],frm[2],50)
        ims.append(im.collections)

    # Create the gif from the plots
    ani = animation.ArtistAnimation(fig, ims, interval=100,repeat=True)

    # Save the final result
    ani.save("diffb.gif", writer='imagemagick')


################################################################################

def gifp(sz):
# Function which takes the 'X' delimited values of the wave, plots and saves
# to a gif. This is not an ideal implementation, more work is required for
# speed up.

    # Initialize
    frames = []
    cur_frame = []
    cur_frame.append([])
    cur_frame.append([])
    cur_frame.append([])

    count = []
    x = []
    y = []
    va = []
    ims = []

    fig = plt.figure()

    file = "time_flow_p.dat"
    x_lab = "X (m)"
    y_lab = "Y (m)"
    title = "Advection-Diffusion \n Concentration of Phage (1/mL)"

    fr = open(file,'r')

    # Read in data
    for line in fr:
        lin = line.split()
        # Save frames
        if lin[0] == "X":
            cur_frame[0] = x
            cur_frame[1] = y
            cur_frame[2] = va
            frames.append(cur_frame)
            cur_frame = []
            cur_frame.append([])
            cur_frame.append([])
            cur_frame.append([])
            x = []
            y = []
            va = []
        else: # Build frames
            x.append(float(lin[0]))
            y.append(float(lin[1]))
            va.append(float(lin[2]))


    # Create the plots
    for k,frm in enumerate(frames):

        plt.title(title)
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)
        im = plt.tricontourf(frm[0],frm[1],frm[2],50)
        ims.append(im.collections)

    # Create the gif from the plots
    ani = animation.ArtistAnimation(fig, ims, interval=100,repeat=True)

    # Save the final result
    ani.save("diffp.gif", writer='imagemagick')


################################################################################

x = []
y = []
C = []

with open("time_mid_p.dat",'r') as f:
    vals = f.read().splitlines()

for line in vals:
    dat = line.split()
    x.append(dat[0])
    y.append(dat[1])
    C.append(dat[2])

for index, val in enumerate(x):
    x[index] = float(val)

for index, val in enumerate(y):
    y[index] = float(val)

for index, val in enumerate(C):
    C[index] = float(val)

plt.tricontourf(x,y,C,30)
plt.colorbar()
plt.title("Concentration Count/mL")
plt.xlabel("cm")
plt.ylabel("cm")
plt.savefig("Finalp.png",bbox_inches='tight')

plt.clf()

x = []
y = []
C = []

with open("time_mid_b.dat",'r') as f:
    vals = f.read().splitlines()

for line in vals:
    dat = line.split()
    x.append(dat[0])
    y.append(dat[1])
    C.append(dat[2])

for index, val in enumerate(x):
    x[index] = float(val)

for index, val in enumerate(y):
    y[index] = float(val)

for index, val in enumerate(C):
    C[index] = float(val)

plt.tricontourf(x,y,C,30)
plt.colorbar()
plt.title("Concentration Count/mL")
plt.xlabel("cm")
plt.ylabel("cm")
plt.savefig("Finalb.png",bbox_inches='tight')



gifb(20)
gifp(20)
