#! /usr/bin/python

# If the input value of 0 is provided to the
# 'timing_comp' binary, the file 'samples.csv'
# is written, which lists planes generated
# in the same way they are for the tests.
# Note, these are not the same planes used,
# but should be statistically similar.
#
# This file plots the normals of the planes on
# the unit sphere.

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

xs = []
ys = []
zs = []
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
with open('samples.csv', newline='\n') as f:
   reader = csv.reader(f,delimiter=',')
   for row in reader:
       normal = [float(row[0]),float(row[1]),float(row[2])]
       xs.append(normal[0])
       ys.append(normal[1])
       zs.append(normal[2])       

ax.scatter(xs, ys, zs, zdir='z', s=20, c=None, depthshade=True)       
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
