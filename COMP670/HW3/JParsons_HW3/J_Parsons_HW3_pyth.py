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

# Python code for HW2, COMP670

# Jon Parsons
# 10-1-19

################################################################################

## Function to handle plotting of the time and space contour plots

def cont_plot(grid):
    # Initialize file list
    file_list = []
    # Define file names, output, axis labeling
    num_file = "num_grid{}.dat".format(grid)
    out_num = "{}pnt_num.png".format(grid)
    title_num = "Wave  ({} Points) \n Numerical".format(grid)
    xlab = "X"
    ylab = "Y"

    file_list.append(num_file)
    # Read in Files

    for item in file_list:
        # Initialize lists
        ln = []
        x = []
        t = []

        # Take in data
        with open(item,'r') as f:
            dat = f.read().splitlines()
        # Split data
        for line in dat:
            lin = line.split()
            t.append(lin[0])
            x.append(lin[1])
            ln.append(lin[2])
        # Convert data to numerical floats
        for index, val in enumerate(x):
            x[index] = float(val)

        for index, val in enumerate(t):
            t[index] = float(val)

        for index, val in enumerate(ln):
            ln[index] = float(val)
        # Plots

        plt.tricontourf(t,x,ln,50)
        clt = plt.colorbar()
        clt.ax.set_title('Amplitude')
        plt.title(title_num)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.savefig(out_num,bbox_inches='tight')


        plt.clf()
################################################################################
def gif(sz):
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

    file = "{}times.dat".format(sz)
    x_lab = "X"
    y_lab = "Y"
    title = "Wave Evolution, {} Points".format(sz)

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
    ani.save("Wave.gif", writer='imagemagick')


################################################################################

# File names known
err_file = "err_pts.dat"

# Number of gridpoints (known)
grid = [20,40,80]

# Error plotting
l2 = [] # holds l2 values
h = [] # holds values of h
nmpts = [] # holds number of gridpoints
pt_err = [] # holds error at a point

################################################################################
# This section handles plotting for the data over time and space
# Outputs are contour plots

#for item in grid:
    #cont_plot(item)

gif(20)
