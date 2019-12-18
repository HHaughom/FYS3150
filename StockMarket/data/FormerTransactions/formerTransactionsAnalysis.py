# Python script for data analysing
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def powerlaw(x, exponent, c, c0):
	return c0 + x**exponent * c

if __name__ == '__main__':

	""" 
	N = 1000
	alpha = 1.0
	lambda = 0.0
	gamma = 1, 2, 3, 4
	"""
	
	alfa1gamma1 = np.loadtxt("alfa1gamma1.dat")
	alfa1gamma2 = np.loadtxt("alfa1gamma2.dat")
	alfa1gamma3 = np.loadtxt("alfa1gamma3.dat")
	alfa1gamma4 = np.loadtxt("alfa1gamma4.dat")
	
	alfa1 = [alfa1gamma1, alfa1gamma2, alfa1gamma3, alfa1gamma4]
	
	""" 
	N = 1000
	alpha = 2.0
	lambda = 0.0
	gamma = 1, 2, 3, 4
	"""
	
	alfa2gamma1 = np.loadtxt("alfa2gamma1.dat")
	alfa2gamma2 = np.loadtxt("alfa2gamma2.dat")
	alfa2gamma3 = np.loadtxt("alfa2gamma3.dat")
	alfa2gamma4 = np.loadtxt("alfa2gamma4.dat")
	
	alfa2 = [alfa2gamma1, alfa2gamma2, alfa2gamma3, alfa2gamma4]
	
	""" 
	N = 1000
	alpha = 1.0
	lambda = 0.5
	gamma = 1, 2, 3, 4
	"""
	
	lambda05alfa1gamma1 = np.loadtxt("lambda05alfa1gamma1.dat")
	lambda05alfa1gamma2 = np.loadtxt("lambda05alfa1gamma2.dat")
	lambda05alfa1gamma3 = np.loadtxt("lambda05alfa1gamma3.dat")
	lambda05alfa1gamma4 = np.loadtxt("lambda05alfa1gamma4.dat")
	
	lambda05alfa1 = [lambda05alfa1gamma1, lambda05alfa1gamma2, lambda05alfa1gamma3, lambda05alfa1gamma4]
	
	""" 
	N = 1000
	alpha = 2.0
	lambda = 0.5
	gamma = 1, 2, 3, 4
	"""
	
	lambda05alfa2gamma1 = np.loadtxt("lambda05alfa2gamma1.dat")
	lambda05alfa2gamma2 = np.loadtxt("lambda05alfa2gamma2.dat")
	lambda05alfa2gamma3 = np.loadtxt("lambda05alfa2gamma3.dat")
	lambda05alfa2gamma4 = np.loadtxt("lambda05alfa2gamma4.dat")
	
	lambda05alfa2 = [lambda05alfa2gamma1, lambda05alfa2gamma2, lambda05alfa2gamma3, lambda05alfa2gamma4]
	
	
	x = np.linspace(0, 20, len(alfa1gamma1))
	labels = ["$\gamma = 1.0$", "$\gamma = 2.0$", "$\gamma = 3.0$", "$\gamma = 4.0$"]
	mark = ['-', '--', '-.', ':']
	
	#Plot distributions 
	
	for para in [alfa1, alfa2, lambda05alfa1, lambda05alfa2]:
		for data in range(len(labels)):
			plt.plot(x, para[data], mark[data], label = labels[data])
			
		plt.yscale("log")
		#plt.xscale("log")
		plt.legend(loc="best", fontsize=13)
		plt.xlabel('m', fontsize=13)
		plt.ylabel('P(m)', fontsize=13)
		plt.show()
		
	
	#Plot tails with power fit 
	
	tailStart = 0.01 #Fraction of the maximum peak where the tail starts
		
	for combination in [alfa1, alfa2, lambda05alfa1, lambda05alfa2]:
		for gamma in range(len(labels)):
			
			data = combination[gamma]
			peak = np.argmax(data)
			tail = np.argmin(np.abs(data[peak]*tailStart - data[peak:])) + peak
			
			plt.scatter(x[tail:], data[tail:], marker='.')
			
			popt = curve_fit(powerlaw, x[tail:], data[tail:], p0 = np.asarray([-10.0,1.0,1.0]))[0]
			popt = curve_fit(powerlaw, x[tail:], data[tail:], p0 = popt)[0]
			plt.plot(x[tail:], powerlaw(x[tail:], *popt), label = (labels[gamma] + ', power: %.2f' % popt[0]))
	
		plt.legend()
		plt.xlabel('m', fontsize=13)
		plt.ylabel('P(m)', fontsize=13)
		plt.show()
	
		
	
