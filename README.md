# HASlib.jl

[![Build Status](https://travis-ci.com/feanor12/HASlib.jl.svg?branch=master)](https://travis-ci.com/feanor12/HASlib.jl)

Library for analysing helium atom scattering data


```julia
# define some constants
mel = 1.6021766208000002e-22 # J
A = 1e-10 # 1 angstrom in meters
a = 4.3*A # lattice constant
Ei = 4.196*mel # helium beam energy
theta = 45.75 # incident angle

# potential parameters
D = 5*mel 
h = 0.1
xi = 0.38/A

# helium 3 mass
m_He3 = 5.0082343984338715e-27 # kg

# beam wave vector
ki = sqrt(Ei/HASlib.h2m(m_He3))

# define a hexagonal grid
grid = HASlib.HexGrid(a)
# set beam direction to 1*a1 0*a2
b1,b2 = HASlib.rotate_reciprocal(grid.l,[1,0])
# get channels
channels = HASlib.get_channels(b1,b2,ki,theta*pi/180,5)
# sort channels and limit to 5 open and 0 closed channels
channels = HASlib.sort_and_cut(channels,5,0)

# calculate scattered intensities
cc(channels,a,h,D,xi,theta,m_He3,(-1*A
```
