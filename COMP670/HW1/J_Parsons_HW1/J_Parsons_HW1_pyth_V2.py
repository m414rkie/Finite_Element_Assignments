#!/usr/bin/python3

from scipy.interpolate import griddata
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# This code takes known output files and plots the data within.

# 1.1 - Replaces the repetitive graphing code with a function



# Python code for HW1, COMP670

# Jon Parsons
# 9-12-19

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
    title_ana = "Wave Evolution over Time ({} Points) \n Analytical".format(grid)
    title_num = "Wave Evolution over Time ({} Points) \n Numerical".format(grid)
    xlab = "Time"
    ylab = "Space"

    file_list.append(num_file)
    file_list.append(ana_file)
    # Read in Files

    flag = 0

    for item in file_list:
        # Initialize lists
        ln = []
        x = []
        t = []

        flag += 1

        with open(item,'r') as f:
            dat = f.read().splitlines()

        for line in dat:
            lin = line.split()
            t.append(lin[0])
            x.append(lin[1])
            ln.append(lin[2])

        for index, val in enumerate(x):
            x[index] = float(val)

        for index, val in enumerate(t):
            t[index] = float(val)

        for index, val in enumerate(ln):
            ln[index] = float(val)

        if flag == 1:
             plt.tricontourf(t,x,ln,30)
             plt.colorbar()
             plt.title(title_num)
             plt.xlabel(xlab)
             plt.ylabel(ylab)
             plt.savefig(out_num,bbox_inches='tight')

             plt.clf()
        elif flag == 2:
             plt.tricontourf(t,x,ln,30)
             plt.colorbar()
             plt.title(title_ana)
             plt.xlabel(xlab)
             plt.ylabel(ylab)
             plt.savefig(out_ana,bbox_inches='tight')

             plt.clf()



################################################################################

# File names known
err_file = "err_grid.dat"

# Number of gridpoints (known)
grid = [10,20,40]

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
    pt_err.append(dat[1])
    h.append(dat[2])
    l2.append(dat[3])

# Reverse order of h. This allows for the plot to be consistent with
# familiar representations of similar data
h.reverse()

# Plot and define values for the point error
plt.xticks(np.arange(0,50, 10))
plt.plot(nmpts, pt_err)
plt.ylabel('Error at (1,2)')
plt.xlabel('Number of Points')
plt.title('Point Error of Numerical Solution')
plt.savefig('err_pnt.png',bbox_inches='tight')

plt.clf()

# Define and plot values for the l2 norm
plt.plot(h, l2)
plt.ylabel('l2 Norm')
plt.xlabel('h')
plt.title('l2 Norm Against h (log-log)')
plt.savefig('l2.png',bbox_inches='tight')

plt.clf()

################################################################################
# This section handles plotting for the data over time and space
# Outputs are contour plots

for item in grid:
    cont_plot(item)
