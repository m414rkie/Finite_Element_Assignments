#!usr/bin/python3
 
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# How does modulo work?

a = int(input("Enter a number: "))
n = int(input("Enter another number: "))

out = a%(n-1)

print("Modulo: ", out)
