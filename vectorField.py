#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inputs
Outputs

Aug. '18
@author: jmobrecht
"""

###############################################################################
# Simulated Field
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
# import pandas as pd
import scipy as sp
# import tensorflow as tf

###############################################################################
# Simulated Field
###############################################################################

# Input Variables
ws  = 6.0
ti  = 0.05
alp = 0.25
th  = 7.5   *   np.pi/180
ph  = 5.5   *   np.pi/180

# Grid Dimensions
nGP = 10    # Number of grid points
nSL = 3     # Number of slices

# Grid Bounds
gMax = nGP + 0.5
gMin = 0.5
ygb, zgb = np.meshgrid(np.linspace(gMin,gMax,2),np.linspace(gMin,gMax,2))

###############################################################################
# Rotor Variables
###############################################################################

# Rotor Dimensions
zHH = nGP/2         # "Hub Height" of the rotor
Rad = 0.9*nGP/2    # Rotor Radius
nr = 7
muNodes = (np.linspace(0,1,nr))**(0.5)
muMidpts = muNodes[1:]-np.diff(muNodes)/2
rv = Rad*muMidpts
nth = 24
thv = np.linspace(360/nth,360,nth)

# Rotor Grid
rR, thR = np.meshgrid(rv,thv)
yR = rR*np.sin(thR*np.pi/180) + (nGP+1)/2
zR = rR*np.cos(thR*np.pi/180) + (nGP+1)/2

# Rotor Perimeter
nP = 100
yP = Rad*np.sin(np.linspace(0,2*np.pi,nP)) + (nGP+1)/2
zP = Rad*np.cos(np.linspace(0,2*np.pi,nP)) + (nGP+1)/2

###############################################################################
# Wind Field
###############################################################################

# Turbulence Field
du = 2 * np.random.random([nSL,nGP,nGP]) - 1
dv = 2 * np.random.random([nSL,nGP,nGP]) - 1
dw = 2 * np.random.random([nSL,nGP,nGP]) - 1

# Spatial Grids
xv = np.linspace(1,nSL,nSL)
yv = np.linspace(1,nGP,nGP)
zv = np.linspace(1,nGP,nGP)
xg, yg, zg = np.meshgrid(xv,yv,zv, indexing='ij')
yg2, zg2 = np.meshgrid(yv,zv)

# Wind Speed
u0 = ws*(zg/zHH)**alp
u = u0*np.cos(th)*np.cos(ph) + ws*ti*du
v = u0*np.sin(th)*np.cos(ph) + ws*ti*dv
w = u0*np.sin(ph)            + ws*ti*dw
s = np.sqrt(u**2 + v**2 + w**2)

###############################################################################
# Calculations
###############################################################################

# Interpolation setup
from scipy.interpolate import interp2d
yRr = yR.ravel()
zRr = zR.ravel()

# Interpolate longitudinal averages of u,v,w,s onto Rotor grid
# Initialize interpolation functions
intkind = 'linear'
fu = interp2d(yg2,zg2,np.mean(u, axis=0),kind=intkind)
fv = interp2d(yg2,zg2,np.mean(v, axis=0),kind=intkind)
fw = interp2d(yg2,zg2,np.mean(w, axis=0),kind=intkind)
fs = interp2d(yg2,zg2,np.mean(s, axis=0),kind=intkind)
# Initialize output arrays
uR = np.zeros((np.size(yRr)))
vR = np.zeros((np.size(yRr)))
wR = np.zeros((np.size(yRr)))
sR = np.zeros((np.size(yRr)))
# Interpolate values onto Rotor grid
for i in range(0,np.size(yRr)):
    uR[i] = fu(yRr[i],zRr[i])
    vR[i] = fv(yRr[i],zRr[i])
    wR[i] = fw(yRr[i],zRr[i])
    sR[i] = fs(yRr[i],zRr[i])

# Calculate wind speed from averaging 'u' values
oMWS = np.mean(uR)

# Calculate horizontal wind speed from 'v' values
oALH = np.arctan2(np.mean(vR),np.mean(uR)) * 180/np.pi
# Calculate horizontal wind speed from 'w' values
oALV = np.arctan2(np.mean(wR),np.mean(uR)) * 180/np.pi

# Calculate shear
from poly1 import poly1
m, b = poly1(np.log(yRr/zHH),np.log(sR))
oMWS2 = np.exp(b)
oSHR = m

###############################################################################
# Plot & Print Stuff
###############################################################################
print()
print('INPUTS:')
print('Horizontal Wind Speed: {:1.4} m/s'.format(ws))
print('Horizontal Alignment: {:1.4} deg'.format(th*180/np.pi))
print('Vertical Alignment: {:1.4} deg'.format(ph*180/np.pi))
print('Wind Shear: {:1.4}'.format(alp))
print('Turbulence Intensity: {:1.4}%'.format(ti*1E2))
print()
print('OUTPUTS:')
print('Horizontal Wind Speed: {:1.4} m/s'.format(oMWS))
print('Horizontal Alignment: {:1.4} deg'.format(oALH))
print('Vertical Alignment: {:1.4} deg'.format(oALV))
print('Wind Shear: {:1.4}'.format(oSHR))
print()

#from mpl_toolkits.mplot3d import axes3d
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.quiver(xg, yg, zg, u, v, w, length=0.5, normalize=True)
#plt.show()

#plt.figure(1)
#plt.plot(yR,zR,'r.',markersize=4)
#plt.plot(yP,zP,'b')
#plt.plot(yg2,zg2,'k.',markersize=1)
#plt.plot([gMin,gMin],[gMin,gMax],'k')
#plt.plot([gMax,gMax],[gMin,gMax],'k')
#plt.plot([gMin,gMax],[gMin,gMin],'k')
#plt.plot([gMin,gMax],[gMax,gMax],'k')
#plt.axis('equal')
#plt.axis('off')

#plt.figure(2)
#plt.plot(yRr,sR,'r.',markersize=4)
#plt.axis('equal')

###############################################################################
# End
###############################################################################