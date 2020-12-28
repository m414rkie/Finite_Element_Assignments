#!/usr/bin/python3

from scipy.interpolate import griddata
import matplotlib as mpl
import matplotlib.pyplot as plt
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
    ana_file = "num_grid{}.dat".format(grid)
    out_ana = "{}pnt_ana.png".format(grid)
    out_num = "{}pnt_num.png".format(grid)
    title_ana = "Wave  ({} Points) \n Analytical".format(grid)
    title_num = "Wave  ({} Points) \n Numerical".format(grid)
    xlab = "X"
    ylab = "Y"

    file_list.append(num_file)
    file_list.append(ana_file)
    # Read in Files

    # Flag determines if plotting numerical or analytical
    flag = 0

    for item in file_list:
        # Initialize lists
        ln = []
        x = []
        t = []

        flag += 1
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
        if flag == 1:
             plt.tricontourf(t,x,ln,50)
             plt.colorbar()
             plt.title(title_num)
             plt.xlabel(xlab)
             plt.ylabel(ylab)
             plt.savefig(out_num,bbox_inches='tight')

             plt.clf()
        elif flag == 2:
             plt.tricontourf(t,x,ln,50)
             plt.colorbar()
             plt.title(title_ana)
             plt.xlabel(xlab)
             plt.ylabel(ylab)
             plt.savefig(out_ana,bbox_inches='tight')

             plt.clf()

################################################################################

# File names known
err_file = "err_pts.dat"

# Number of gridpoints (known)
grid = [15,30,60]

# Error plotting
l2 = [] # holds l2 values
h = [] # holds values of h
nmpts = [] # holds number of gridpoints
pt_err = [] # holds error at a point

# Read in the error file
with open(err_file,'r') as f:
    errdat = f.read().splitlines()

# Convert data to usable form
for line in errdat:
    dat = line.split()
    nmpts.append(float(dat[0]))
    l2.append((float(dat[1])))

# familiar representations of similar data
l2_pl = [] # For Plotting
l2_ln = [] # For poly fit
h_ln = [] # For poly fit
h_sq = [] # For comparison

for index, val in enumerate(l2):
    l2_pl.append((np.float(val)))
    l2_ln.append(math.log(np.float(val)))

for i, pt in enumerate(nmpts):
    h.append((np.float(pt)))
    h_sq.append(np.float(pt)**2)
    h_ln.append(math.log(np.float(pt)))

fit = np.polyfit(h_ln,l2_ln,1)

print("Polynomial fit of error on a natural log scale")
print("y = {}ln(err) + {}".format(fit[0],fit[1]))

# Define and plot values for the l2 norm
plt.loglog(h, l2_pl,label="y = {}ln(err) {}".format(fit[0],fit[1]))
plt.loglog(h, l2_pl,'bo')
plt.loglog(h,h_sq,'r--',label='h^2 ')
plt.legend(loc='upper left')
plt.ylabel('l2 Norm')
plt.xlabel('h')
plt.title('l2 Norm Against h (Ln - Ln)')
plt.savefig('l2ll.png',bbox_inches='tight')

plt.clf()

################################################################################
# This section handles plotting for the data over time and space
# Outputs are contour plots

for item in grid:
    cont_plot(item)
