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
fp=open("ref_sol.dat","w")
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
fp=open("samp_pts.dat","w")
for i in xrange(ndp):
    w=np.dot(v,x_meas[i,:])
    fp.write(" ".join([str(a) for a in [time[(i+1)*q],x_meas[i,0],x_meas[i,1],x_meas[i,2],w[0],w[1],w[2]]])+"\n")
fp.close()

# Define a factor controlling the increase in variance in the model (i.e. the
# term Q in the Kalman notation). This is based on an estimate of numerical
# integration error. Python's 'odeint' uses absolute/relative error tolerances
# of ~1e-8 by default. The factor below is based on some padding, plus
# solutions with a scale on the order of ~10.
ifac=ndp*1e-6/tmax

# Kalman filter ODE system: mean movement plus variance update
def dkalman(x,t):

    # Construct the symmetric Kalman P variance matrix using terms 3 to 8
    P=np.array([[x[3],x[4],x[5]],[x[4],x[6],x[7]],[x[5],x[7],x[8]]])

    # Compute dP/dt using the Jacobian of the Chua model, plus the general Q
    # term from the inherent model variance.
    F=np.array([[alpha*df(x[0]),alpha,0],[1,-1,1],[0,-beta,0]])
    S=np.dot(F,P)+np.dot(P,F.T)+ifac*ifac*np.eye(3)

    # Return a vector of updates to the mean and to the variance
    return np.concatenate((deriv(x,t),[S[0,0],S[0,1],S[0,2],S[1,1],S[1,2],S[2,2]]))

# Initialize the mean and variance of the solution
y=np.empty((n,9))
y[0,:]=np.array([0.3,0.2,0.3,sig2,0,0,sig2,0,sig2])
y[0,0:3]+=np.random.normal(0.,sig,(3))

# Loop over the available data points
for i in xrange(ndp):

    # Integrate the Kalman ODE system (mean & variance) up to the next data
    # point
    l=i*q;m=l+q
    y[l:m+1,:]=odeint(dkalman,y[l,:],time[l:m+1])

    # Construct the symmetric Kalman P variance matrix using terms 3 to 8
    P=np.array([[y[m,3],y[m,4],y[m,5]],[y[m,4],y[m,6],y[m,7]],[y[m,5],y[m,7],y[m,8]]])

    # Do the Kalman updates to the mean and variance based on the data
    K=np.dot(P,np.linalg.inv(P+sig2*np.eye(3)))
    y[m,0:3]+=np.dot(K,x_meas[i,:]-y[m,0:3])
    P-=np.dot(K,P)

    # Update terms 3 to 8 using the new P matrix
    y[m,3]=P[0,0];y[m,4]=P[0,1];y[m,5]=P[0,2]
    y[m,6]=P[1,1];y[m,7]=P[1,2];y[m,8]=P[2,2]

# Save the mean and variance of the Kalman filter solution
fp=open("kal_sol.dat","w")
for i in xrange(n):
    w=np.dot(v,y[i,0:3])
    fp.write(" ".join([str(a) for a in np.concatenate(([time[i]],y[i,:],w))])+"\n")
fp.close()
