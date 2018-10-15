#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import *

# Seed the random number generator, so that the code will return the same
# output when run twice. This can help with seeing how the results change to
# minor parameter tweaks.
np.random.seed(1)

# Read in axes that were previously stored by chua.py
v=np.load("axes.npy")

# Constants in model
alpha=10
beta=14

# Nonlinear function in model, and its derivative
def f(u):
    return u*(1/6.-0.0625*u*u)
def df(u):
    return 1/6.-0.1875*u*u

# Function describing the Chua chaotic ODE system
def deriv(x,t):
    return np.array([alpha*(x[1]+f(x[0])),x[0]-x[1]+x[2],-beta*x[1]])

# Calculate reference solution
tmax=200.
n=20001
time=np.linspace(0,tmax,n)
xinit=np.array([0.3,0.2,0.3])
x=odeint(deriv,xinit,time)

# Save reference solution
fp=open("ref_sol2.dat","w")
for i in xrange(n):
    w=np.dot(v,x[i,:])
    fp.write(" ".join([str(a) for a in [time[i],x[i,0],x[i,1],x[i,2],w[0],w[1],w[2]]])+"\n")
fp.close()

# Compute data points by sampling the reference solution at intervals of size
# q
q=50
ndp=(n-1)/q
x_meas=np.empty((ndp,3))
for i in range(ndp):
    x_meas[i,:]=x[(i+1)*q,:]

# Add Gaussian noise with standard deviation sig
sig=7e-2
sig2=sig*sig
x_meas+=np.random.normal(0.,sig,(ndp,3))

# Save data points
fp=open("samp_pts2.dat","w")
for i in xrange(ndp):
    w=np.dot(v,x_meas[i,:])
    fp.write(" ".join([str(a) for a in [time[(i+1)*q],x_meas[i,0],x_meas[i,1],x_meas[i,2],w[0],w[1],w[2]]])+"\n")
fp.close()

# Define a factor controlling the increase in variance in the model (i.e. the
# term Q in the Kalman notation). This is based on an estimate of numerical
# integration error. Python's 'odeint' uses absolute/relative error tolerances
# of ~1e-8 by default. The factor below is based on some padding, plus
# solutions with a scale on the order of ~10.
# Kalman filter ODE system: mean movement plus variance update
ifac=ndp*1e-5/tmax

# Kalman filter: mean movement plus variance update
def dkalman(x,t):

    # Construct the symmetric Kalman P variance matrix using terms 5 to 24.
    # This makes just uses the reshape command and has some redundancy because
    # it solves for all 5x5 components and doesn't exploit symmetry.
    P=x[5:].reshape((5,5))

    # Compute dP/dt using the Jacobian of the Chua model, plus the general Q
    # term from the inherent model variance.
    F=np.array([[x[3]*df(x[0]),x[3],0,x[1]+f(x[0]),0],[1,-1,1,0,0],[0,-x[4],0,0,-x[1]],[0,0,0,0,0],[0,0,0,0,0]])
    S=np.dot(F,P)+np.dot(P,F.T)+np.diag([ifac*ifac,ifac*ifac,ifac*ifac,ifac*ifac,ifac*ifac])

    # Return a vector of updates to the mean and to the variance
    return np.concatenate(([x[3]*(x[1]+f(x[0])),x[0]-x[1]+x[2],-x[4]*x[1],0,0],S.reshape((25))))

# Set up for Kalman filter: initial condition, plus the six independent
# components of the Kalman filter
# Initialize the mean and variance of the solution
y=np.empty((n,30))
y[0,0:5]=np.array([0.3,0.2,0.3,9.6,13.6])
y[0,0:3]+=np.random.normal(0.,sig,(3))
y[0,5:]=np.diag([sig*sig,sig*sig,sig*sig,1.,1.]).reshape((25))

# Define the H matrix describing what measurements are made. Here, only the
# three ODE components are viewed. The model parameters aren't visible
H=np.array([[1.,0.,0.,0.,0.],[0.,1.,0.,0.,0.],[0.,0.,1.,0.,0.]])

# Loop over the available data points
for i in xrange(ndp):

    # Integrate the Kalman ODE system (mean & variance) up to the next data
    # point
    l=i*q;m=l+q
    y[l:m+1,:]=odeint(dkalman,y[l,:],time[l:m+1])

    # Do the Kalman updates to the mean and variance based on the data
    P=y[m,5:].reshape(5,5)
    K=np.dot(np.dot(P,H.T),np.linalg.inv(P[0:3,0:3]+sig*sig*np.eye(3)))
    y[m,0:5]+=np.dot(K,x_meas[i,:]-y[m,0:3])
    P-=np.dot(np.dot(K,H),P)
    y[m,5:]=P.reshape((25))

# Save the mean and variance of the Kalman filter solution
fp=open("kal_sol2.dat","w")
for i in xrange(n):
    w=np.dot(v,y[i,0:3])
    fp.write(" ".join([str(a) for a in np.concatenate(([time[i]],y[i,:],w))])+"\n")
fp.close()
