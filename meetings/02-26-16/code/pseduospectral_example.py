#code adapted from http://www.scientificpython.net/
import numpy as np
import scipy as sp
from operator import mul
from scipy.linalg import solve
import matplotlib.pyplot as plt


def diffmat(x): # x is an ordered array of grid points
	n = sp.size(x)
	e = sp.ones((n,1))
	Xdiff = sp.outer(x,e)-sp.outer(e,x)+sp.identity(n)
	xprod = -reduce(mul,Xdiff) # product of rows
	W = sp.outer(1/xprod,e)
	D = W/sp.multiply(W.T,Xdiff)
	d = 1-sum(D)
	for k in range(0,n):  # Set diagonal elements
		D[k,k] = d[k]
	return -D.T

def true_solution(x):
	y_sol = np.zeros(len(x))
	return np.array([x[i]*np.e**(-x[i]) for i in xrange(len(x))])

if __name__=='__main__':

	n = 20 #Try increasing this
	x = 0.5*(1-np.cos(np.pi*np.arange(n+1)/(n)))
	D = diffmat(x)
	print 'Condition number of D is: ',np.linalg.cond(D)
	I = np.identity(n+1)
	A = D + I
	print A
	plt.imshow(A,interpolation='none')
	plt.show()
	print 'Condition number of A is: ',np.linalg.cond(A)

	f = np.exp(-x)
	y = np.zeros(n+1)
	y[1:] = solve(A[1:,1:],f[1:])
	print y
	true_x = np.linspace(0,1,100)
	true_y = true_solution(true_x)
	true_y2 = true_solution(x)
	plt.plot(true_x,true_y,'b-',label='True Solution')
	plt.plot(x,y,'ro')
	plt.plot(x,y,'k--',label='Collocation Method')
	plt.legend(loc=2)
	plt.show()

	plt.plot(x,y-true_y2,'ko',label='Error')
	plt.legend(loc=2)
	plt.show()

	##Lets investigate condition number
	ns = range(3,25)
	conds = []
	for n in ns:
		x = 0.5*(1-np.cos(np.pi*np.arange(n+1)/(n)))
		D = diffmat(x)
		I = np.identity(n+1)
		A = D + I
		print n,np.linalg.cond(A)
		conds.append(np.linalg.cond(A))
	plt.plot(ns,conds,'ko')
	plt.yscale('log')
	plt.show()





