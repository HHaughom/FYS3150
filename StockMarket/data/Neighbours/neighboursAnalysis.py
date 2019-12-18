# Python script for data analysing
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def powerlaw(x, exponent, c, c0):
	return c0 + x**exponent * c

if __name__ == '__main__':

	""" 
	#N = 1000
	#lambda = 0
	#alpha = 0, 0.5, 1.0, 1.5, 2.0 
	"""
	
	N1000alfa0lambda0 = np.loadtxt("N1000alfa0lambda0.dat")
	N1000alfa05lambda0 = np.loadtxt("N1000alfa05lambda0.dat")
	N1000alfa1lambda0 = np.loadtxt("N1000alfa1lambda0.dat")
	N1000alfa15lambda0 = np.loadtxt("N1000alfa15lambda0.dat")
	N1000alfa2lambda0 = np.loadtxt("N1000alfa2lambda0.dat")
	
	N1000lambda0 = [N1000alfa05lambda0, N1000alfa1lambda0, N1000alfa15lambda0, N1000alfa2lambda0]
	
	""" 
	#N = 1000
	#lambda = 0.5
	#alpha = 0.5, 1.0, 1.5, 2.0
	"""
	
	N1000alfa05lambda05 = np.loadtxt("N1000alfa05lambda05.dat")
	N1000alfa1lambda05 = np.loadtxt("N1000alfa1lambda05.dat")
	N1000alfa15lambda05 = np.loadtxt("N1000alfa15lambda05.dat")
	N1000alfa2lambda05 = np.loadtxt("N1000alfa2lambda05.dat")
	
	N1000lambda05 = [N1000alfa05lambda05, N1000alfa1lambda05, N1000alfa15lambda05, N1000alfa2lambda05]
	
	""" 
	#N = 500
	#lambda = 0
	#alpha = 0.5, 1.0, 1.5, 2.0
	"""
	
	N500alfa05lambda0 = np.loadtxt("N500alfa05lambda0.dat")
	N500alfa1lambda0 = np.loadtxt("N500alfa1lambda0.dat")
	N500alfa15lambda0 = np.loadtxt("N500alfa15lambda0.dat")
	N500alfa2lambda0 = np.loadtxt("N500alfa2lambda0.dat")
	
	N500lambda0 = [N500alfa05lambda0, N500alfa1lambda0, N500alfa15lambda0, N500alfa2lambda0]
	
	""" 
	#N = 500
	#lambda = 0.5
	#alpha = 0.5, 1.0, 1.5, 2.0
	"""
	
	N500alfa05lambda05 = np.loadtxt("N500alfa05lambda05.dat")
	N500alfa1lambda05 = np.loadtxt("N500alfa1lambda05.dat")
	N500alfa15lambda05 = np.loadtxt("N500alfa15lambda05.dat")
	N500alfa2lambda05 = np.loadtxt("N500alfa2lambda05.dat")
	
	N500lambda05 = [N500alfa05lambda05, N500alfa1lambda05, N500alfa15lambda05, N500alfa2lambda05]
	
	x = np.linspace(0, 20, len(N1000alfa0lambda0))
	lab = [r'$\alpha = 0.5$', r'$\alpha = 1.0$', r'$\alpha = 1.5$', r'$\alpha = 2.0$']
	mark = ['-', '--', '-.', ':']
	
	#Plot distributions 
	
	i = 1
	for para in [N1000lambda0, N1000lambda05, N500lambda0, N500lambda05]:
		plt.figure(i)
		
		for data in range(len(lab)):
			plt.plot(x, para[data], mark[data], label = lab[data])
		
		#plt.yscale("log")
		#plt.xscale("log")
		plt.legend(loc="best", fontsize=13)
		plt.xlabel('m', fontsize=13)
		plt.ylabel('P(m)', fontsize=13)
		plt.show()
		i += 1
		
		
	#Plot tails with powerfit
	
	tailStart = 0.01 #Fraction of the maximum peak where the tail starts
	
	for combination in [N1000lambda0, N1000lambda05, N500lambda0, N500lambda05]:
		plt.figure(i)
		for alpha in range(len(lab)):
			
			data = combination[alpha]
			peak = np.argmax(data)
			tail = np.argmin(np.abs(data[peak]*tailStart - data[peak:])) + peak
			
			plt.scatter(x[tail:], data[tail:], marker='.')
			
			popt = curve_fit(powerlaw, x[tail:], data[tail:], p0 = np.asarray([-5,1.0,1.0]))[0]
			popt = curve_fit(powerlaw, x[tail:], data[tail:], p0 = popt)[0]
			plt.plot(x[tail:], powerlaw(x[tail:], *popt), label = (lab[alpha] + ', power: %.2f' % popt[0]))
	
		plt.legend()
		plt.xlabel('m', fontsize=13)
		plt.ylabel('P(m)', fontsize=13)
		plt.show()
			