# Python script for data analysing

import math as math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def w(m, m0):
	"""Gibbs distribution"""
	beta = 1.0/m0
	return beta*np.exp(-beta*m)

def n(saving):
    return 1 + 3*saving/(1-saving)

def a(n):
    return n**n/math.gamma(n)

def Pn(x, saving):
	"""Gamma distribution"""
	var = n(saving)
	return a(var)*x**(var-1)*np.exp(-var*x)
	
def powerlaw(x, exponent, c, c0):
	return c0 + x**exponent * c

if __name__ == '__main__':

	#Loads distribution data
	lambda0 = np.loadtxt("lambda0.dat")
	lambda025 = np.loadtxt("lambda025.dat")
	lambda05 = np.loadtxt("lambda05.dat")
	lambda09 = np.loadtxt("lambda09.dat")
	x = np.linspace(0,10,len(lambda0))

	#Plots raw data points
	plt.figure(1)
	plt.scatter(x, lambda0, marker=".", label="$\lambda = 0$")
	plt.scatter(x, lambda025, marker="+", label="$\lambda = 0.25$")
	plt.scatter(x, lambda05, marker="x", label="$\lambda = 0.5$")
	plt.scatter(x, lambda09, marker="1", label="$\lambda = 0.9$")
	
	#Calculate analytical distributions
	y1 = w(x,1.0)
	y2 = Pn(x, 0.25)
	y3 = Pn(x, 0.5)
	y4 = Pn(x, 0.9)
	
	#Normalise analytical distributions
	y1 /= np.sum(y1)
	y2 /= np.sum(y2)
	y3 /= np.sum(y3)
	y4 /= np.sum(y4)
	
	#Plot analytical distributions
	plt.plot(x,y1)
	plt.plot(x,y2)
	plt.plot(x,y3)
	plt.plot(x,y4)
	plt.yscale("log")
	plt.xscale("log")
	plt.legend(loc="best")
	plt.xlabel('m', fontsize=13)
	plt.ylabel('P(m)', fontsize=13)
	plt.show()
	
	#Plot tails with powerfit
	plt.figure(2)
	
	lab = ['$\lambda = 0.0$', '$\lambda = 0.25$', '$\lambda = 0.5$', '$\lambda = 0.9$']

	tailStart = 0.01 #Fraction of the maximum peak where the tail starts
	j = 0
	for lambda_ in [lambda0, lambda025, lambda05, lambda09]:
		peak = np.argmax(lambda_)
		tail = np.argmin(np.abs(lambda_[peak]*tailStart - lambda_[peak:])) + peak
		
		plt.scatter(x[tail:], lambda_[tail:], marker = '.')
		
		popt = curve_fit(powerlaw, x[tail:], lambda_[tail:], p0 = np.asarray([-4,1.0,0.05]))[0]
		popt = curve_fit(powerlaw, x[tail:], lambda_[tail:], p0 = popt)[0]
		plt.plot(x[tail:], powerlaw(x[tail:], *popt), label = (lab[j] + ', power: %.2f' % popt[0]))
		
		j += 1
		
	plt.legend()
	plt.xlabel('m', fontsize=13)
	plt.ylabel('P(m)', fontsize=13)
	plt.show()