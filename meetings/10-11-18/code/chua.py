#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import *

# Constants in model
alpha=10.4
beta=14

# Nonlinear function in model
def f(v):
    return v*(1/6.-0.0625*v*v)

# Function describing the Chua chaotic ODE system
def deriv(x,t):
    return np.array([alpha*(x[1]+f(x[0])),x[0]-x[1]+x[2],-beta*x[1]])

# Solve ODE using the "odeint" library in SciPy
tmax=250.
n=20001
time=np.linspace(0,tmax,n)

# Initial conditions
xinit=np.array([0.3,0.2,0.3])
x=odeint(deriv,xinit,time)

# Perform SVD to find a better set of axes
(u,s,v)=np.linalg.svd(x,full_matrices=False)

# Scale the axes by the singular values. Save for use with other programs.
fac=sqrt(n)
v=np.dot(fac*np.diag([1./s[0],1./s[1],1./s[2]]),v)
np.save("axes.npy",v)

# Save the results
fp=open("chua_sol.dat","w")
for i in xrange(n):
    w=np.dot(v,x[i,:])
    fp.write(" ".join([str(a) for a in time[i],x[i,0],x[i,1],x[i,2],w[0],w[1],w[2]]+"\n"))
fp.close()
