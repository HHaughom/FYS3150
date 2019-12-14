import numpy as np
import matplotlib.pyplot as plt
import sys
import collections
from scipy.signal import savgol_filter

def readExpectationValues(filname):
	data = open(filname).readlines()[1:]
	N = len(data)
	
	mcs = np.zeros(N)
	E = np.zeros(N)
	Mabs = np.zeros(N)
	Cv = np.zeros(N)
	chiAbs = np.zeros(N)
	
	i = 0
	
	for line in data:
		mcs[i] = float(line.split()[0])
		E[i] = float(line.split()[1])
		Mabs[i] = float(line.split()[2]) 
		Cv[i] = float(line.split()[3]) 
		chiAbs[i] = float(line.split()[4])
			
		i += 1
	
	return mcs, E, Mabs, Cv, chiAbs
		
def readAcceptRate(file):
	data = open(file).readlines()[1:]
	N = len(data)
	
	mcs = np.zeros(N)
	accepted = np.zeros(N)
	
	i = 0
	
	for line in data:
		mcs[i] = float(line.split()[0])
		accepted[i] = float(line.split()[1])
			
		i += 1
	
	return mcs, accepted
		
def readEnergyData(file):

	data = open(file).readlines()[1:]
	N = len(data)-1

	E = [float(line.split()[0]) for line in data]
	Ecount = collections.Counter(E)
	keys = Ecount.keys()

	n = len(keys)
	
	E = np.zeros(n)
	frequency = np.zeros_like(E)

	i = 0
	for key in keys:
		E[i] = key
		frequency[i] = Ecount[key]
		i += 1

	Ecopy = E.copy()
	fcopy = frequency.copy()

	Esorted = np.zeros_like(E)
	fsorted = np.zeros_like(E)

	for j in range(n):
		iMin = Ecopy.argmin()
		Esorted[j] = Ecopy[iMin]
		fsorted[j] = fcopy[iMin]
		
		Ecopy = np.delete(Ecopy, iMin)
		fcopy = np.delete(fcopy, iMin)
		
	return Esorted, fsorted, N

def readExpectationTemperatureValues(file):
	data = open(file).readlines()[2:]
	N = len(data)
	
	T = np.zeros(N)
	E = np.zeros(N)
	Mabs = np.zeros(N)
	Cv = np.zeros(N)
	chiAbs = np.zeros(N)
	
	i = 0
	
	for line in data:
		T[i] = float(line.split()[0])
		E[i] = float(line.split()[1])
		Mabs[i] = float(line.split()[2]) 
		Cv[i] = float(line.split()[3]) 
		chiAbs[i] = float(line.split()[4])
			
		i += 1
	
	return T, E, Mabs, Cv, chiAbs

def smooth(data, windowSize = 23,  degree = 5):
	smoothed = savgol_filter(data, windowSize, degree) # window size 51, polynomial order 5
	return smoothed

if sys.argv[1] == "EnergyMC":

	#Plot Energy over MC cycles, task 4c)

	mcs, E24r, Mabs24r, Cv24r, ChiAbs24r = readExpectationValues('L20T24random.txt')
	E24o, Mabs24o, Cv24o, ChiAbs24o = readExpectationValues('L20T24ordered.txt')[1:]
	E1r, Mabs1r, Cv1r, ChiAbs1r = readExpectationValues('L20T1random.txt')[1:]
	E1o, Mabs1o, Cv1o, ChiAbs1o = readExpectationValues('L20T1ordered.txt')[1:]

	fig = plt.figure(figsize=(10.5,3))
	ax = fig.add_subplot(121)
	ax.grid(True)
	plt.title("T=2.4")
	plt.plot(mcs, E24r, linewidth = 3, label = 'Random')
	plt.plot(mcs, E24o, linewidth = 3, label = 'Ordered')
	plt.xscale('log')
	ax.set_xlabel('Monte Carlo cycles', fontsize=12)
	ax.set_ylabel('<E>', fontsize=12)	
	ax.legend(loc=1, prop={'size': 13})

	ax = fig.add_subplot(122)
	ax.grid(True)
	plt.title("T=1.0")
	plt.plot(mcs, E1r, linewidth = 3, label = 'Random')
	plt.plot(mcs, E1o, linewidth = 3, label = 'Ordered')
	plt.xscale('log')
	ax.set_xlabel('Monte Carlo cycles', fontsize=12)
	ax.set_ylabel('<E>', fontsize=12)	
	ax.legend(loc=1, prop={'size': 13})
	fig.savefig('E.png', bbox_inches="tight")

if sys.argv[1] == "MagnetismMC":
	
	#Plot Magnetism over MC cycles, task 4c)
	
	mcs, E24r, M24r, Mabs24r, Cv24r, Chi24r, ChiAbs24r = readExpectationValues('L20T24random.txt')
	E24o, M24o, Mabs24o, Cv24o, Chi24o, ChiAbs24o = readExpectationValues('L20T24ordered.txt')[1:]
	E1r, M1r, Mabs1r, Cv1r, Chi1r, ChiAbs1r = readExpectationValues('L20T1random.txt')[1:]
	E1o, M1o, Mabs1o, Cv1o, Chi1o, ChiAbs1o = readExpectationValues('L20T1ordered.txt')[1:]


	fig = plt.figure(figsize=(10.5,3))
	ax = fig.add_subplot(121)
	ax.grid(True)
	plt.title("T=2.4")
	plt.plot(mcs, Mabs24r, linewidth = 3, label = 'Random')
	plt.plot(mcs, Mabs24o, linewidth = 3, label = 'Ordered')
	plt.yscale('linear')
	plt.xscale('log')
	ax.set_xlabel('Monte Carlo cycles', fontsize=12)
	ax.set_ylabel('<|M|>', fontsize=12)	
	ax.legend(loc=4, prop={'size': 13})

	ax = fig.add_subplot(122)
	ax.grid(True)
	plt.title("T=1.0")
	plt.plot(mcs, Mabs1r, linewidth = 3, label = 'Random')
	plt.plot(mcs, Mabs1o, linewidth = 3, label = 'Ordered')
	plt.yscale('linear')
	plt.xscale('log')
	ax.set_xlabel('Monte Carlo cycles', fontsize=12)
	ax.set_ylabel('<|M|>', fontsize=12)	
	ax.legend(loc=4, prop={'size': 13})
	fig.savefig('M.png', bbox_inches="tight")
	fig = plt.figure(figsize=(4,3))
	
if sys.argv[1] == "AR":

	#Plot state acceptance rate in Metropolis over MC cycles, 4c)

	mcs, accepted1 = readAcceptRate("ET1accepted.txt")
	mcs, accepted2 = readAcceptRate("ET24accepted.txt")

	fig = plt.figure(figsize=(4,3))
	ax = fig.add_subplot(111)
	ax.grid(True)
	plt.plot(mcs, accepted1/mcs, linewidth = 3, label = 'T = 1.0')
	plt.plot(mcs, accepted2/mcs, linewidth = 3, label = 'T = 2.4')
	plt.yscale('linear')
	plt.xscale('log')
	ax.set_xlabel('Monte Carlo cycles', fontsize=12)
	ax.set_ylabel('Acceptance rate', fontsize=12)
	ax.legend(loc='best', prop={'size': 13})
	fig.savefig('acceptanceRate.png', bbox_inches="tight")

if sys.argv[1] == "Evariance":

	#Plot Energy variance, task 4d)

	mcs, E24r, M24r, Mabs24r, Cv24r, Chi24r, ChiAbs24r = readExpectationValues('L20T24random.txt')
	E1r, M1r, Mabs1r, Cv1r, Chi1r, ChiAbs1r = readExpectationValues('L20T1random.txt')[1:]

	fig = plt.figure(figsize=(4,3))
	ax = fig.add_subplot(111)
	ax.grid(True)
	
	plt.plot(mcs, Cv1r/mcs, linewidth = 3, label = "T = 1.0")
	plt.plot(mcs, Cv24r*2.4**2/mcs, linewidth = 3, label = "T = 2.4")

	plt.yscale('log')
	plt.xscale('log')
	
	ax.set_xlabel('Monte Carlo cycles', fontsize=12)
	ax.set_ylabel('var(E)', fontsize=12)	

	ax.legend(loc='best', prop={'size': 13})
	fig.savefig('Evariance.png', bbox_inches="tight")
	
if sys.argv[1] == "PE":

	#Plot distribution of accepted energies

	E1, f1, N = readEnergyData("ET1data.txt")
	E2, f2 = readEnergyData("ET24data.txt")[:-1]

	fig = plt.figure(figsize=(12,3))
	ax = fig.add_subplot(121)
	ax.grid(True)
	plt.title("T=1.0")
	plt.plot(E1, f1/N, linewidth = 3)
	ax.set_xlabel('E', fontsize=12)
	ax.set_ylabel('P(E)', fontsize=12)

	ax = fig.add_subplot(122)
	ax.grid(True)
	plt.title("T=2.4")
	plt.plot(E2, f2/N, linewidth = 3)
	ax.set_xlabel('E', fontsize=12)
	ax.set_ylabel('P(E)', fontsize=12)
	fig.savefig('PE.png', bbox_inches="tight")

if sys.argv[1] == "ExpT":
	
	#plot expectation values over temperature, task e)
	
	file20 = "L20lang.txt"
	file40 = "L40lang.txt"
	file60 = "L60lang.txt"
	file80 = "L80lang.txt"
	file100 = "L100lang.txt" 

	T, EL20, MabsL20, CvL20,  chiAbsL20 = readExpectationTemperatureValues(file20)	
	T, EL40,  MabsL40, CvL40,  chiAbsL40 = readExpectationTemperatureValues(file40)
	EL60, MabsL60, CvL60, chiAbsL60 = readExpectationTemperatureValues(file60)[1:]
	EL80, MabsL80, CvL80, chiAbsL80 = readExpectationTemperatureValues(file80)[1:]
	EL100, MabsL100, CvL100, chiAbsL100 = readExpectationTemperatureValues(file100)[1:]
		
	fig = plt.figure(figsize=(4,3))
	ax = fig.add_subplot(111)
	ax.grid(True)
	plt.plot(T, EL20 , linewidth = 3, label = 'L=20')
	#plt.plot(T, smooth(EL20), linewidth = 3, label = 'L=20')
	plt.plot(T, EL40 , linewidth = 3, label = 'L=40')
	#plt.plot(T, smooth(EL40), linewidth = 3, label = 'L=40')
	plt.plot(T, EL60 , linewidth = 3, label = 'L=60')
	#plt.plot(T, smooth(EL60), linewidth = 3, label = 'L=60')
	plt.plot(T, EL80 , linewidth = 3, label = 'L=80')
	#plt.plot(T, smooth(EL80), linewidth = 3, label = 'L=80')
	plt.plot(T, EL100 , linewidth = 3, label = 'L=100')
	#plt.plot(T, smooth(EL100), linewidth = 3, label = 'L=100')

	ax.set_xlabel('T', fontsize=12)
	ax.set_ylabel('<E>', fontsize=12)
	ax.legend(loc='best', prop={'size': 13})
	fig.savefig('E_T.png', bbox_inches="tight")

	fig = plt.figure(figsize=(4,3))
	ax = fig.add_subplot(111)
	ax.grid(True)
	#plt.plot(T, MabsL20 , linewidth = 3, label = 'L=20')
	plt.plot(T, smooth(MabsL20), linewidth = 3, label = 'L=20')
	#plt.plot(T, MabsL40 , linewidth = 3, label = 'L=40')
	plt.plot(T, smooth(MabsL40), linewidth = 3, label = 'L=40')
	#plt.plot(T, MabsL60 , linewidth = 3, label = 'L=60')
	plt.plot(T, smooth(MabsL60), linewidth = 3, label = 'L=60')
	#plt.plot(T, MabsL80 , linewidth = 3, label = 'L=80')
	plt.plot(T, smooth(MabsL80), linewidth = 3, label = 'L=80')
	#plt.plot(T, MabsL100 , linewidth = 3, label = 'L=100')
	plt.plot(T, smooth(MabsL100), linewidth = 3, label = 'L=100')

	ax.set_xlabel('T', fontsize=12)
	ax.set_ylabel('<|M|>', fontsize=12)
	ax.legend(loc='best', prop={'size': 11})
	fig.savefig('Mabs_T.png', bbox_inches="tight")

	fig = plt.figure(figsize=(4,3))
	ax = fig.add_subplot(111)
	ax.grid(True)
	#plt.plot(T, CvL20, linewidth = 3, label = 'L=20')
	plt.plot(T, smooth(CvL20), linewidth = 3, label = 'L=20')
	#plt.plot(T, CvL40 , linewidth = 3, label = 'L=40')
	plt.plot(T, smooth(CvL40), linewidth = 3, label = 'L=40')
	#plt.plot(T, CvL60 , linewidth = 3, label = 'L=60')
	plt.plot(T, smooth(CvL60), linewidth = 3, label = 'L=60')
	#plt.plot(T, CvL80 , linewidth = 3, label = 'L=80')
	plt.plot(T, smooth(CvL80), linewidth = 3, label = 'L=80')
	#plt.plot(T, CvL100 , linewidth = 3, label = 'L=100')
	plt.plot(T, smooth(CvL100), linewidth = 3, label = 'L=100')
	ax.set_xlabel('T', fontsize=12)
	ax.set_ylabel('Cv', fontsize=12)
	ax.legend(loc='best', prop={'size': 8})
	fig.savefig('Cv_T.png', bbox_inches="tight")

	fig = plt.figure(figsize=(4,3))
	ax = fig.add_subplot(111)
	ax.grid(True)
	
	#plt.plot(T, chiAbsL20 , linewidth = 3, label = 'L=20')
	plt.plot(T, smooth(chiAbsL20), linewidth = 3, label = 'L=20')
	#plt.plot(T, chiAbsL40 , linewidth = 3, label = 'L=40')
	plt.plot(T, smooth(chiAbsL40), linewidth = 3, label = 'L=40')
	#plt.plot(T, chiAbsL60 , linewidth = 3, label = 'L=60')
	plt.plot(T, smooth(chiAbsL60), linewidth = 3, label = 'L=60')
	#plt.plot(T, chiAbsL80 , linewidth = 3, label = 'L=80')
	plt.plot(T, smooth(chiAbsL80), linewidth = 3, label = 'L=80')
	#plt.plot(T, chiAbsL100 , linewidth = 3, label = 'L=100')
	plt.plot(T, smooth(chiAbsL100), linewidth = 3, label = 'L=100')
	ax.set_xlabel('T', fontsize=12)
	ax.set_ylabel('Chi', fontsize=12)
	ax.legend(loc='best', prop={'size': 11})
	fig.savefig('ChiAbs_T.png', bbox_inches="tight")

if sys.argv[1] == "Tinfty":

	#estimate from peaks in graphs in e)

	L = np.asarray([40, 60, 80, 100])
	Tc = np.asarray([2.30 , 2.315, 2.36, 2.37])
	
	TcInfty, b = np.polyfit(L, Tc*L, 1)

	fig = plt.figure(figsize=(4,3))
	ax = fig.add_subplot(111)
	ax.grid(True)
	ax.set_xlabel('L', fontsize=12)
	ax.set_ylabel('L*Tc', fontsize=12)
	plt.plot(L, TcInfty*L+b, linewidth =3, label = 'Best-fit')
	plt.plot(L, Tc*L, 'o', color = 'orange', label = 'Data')
	ax.legend(loc='best', prop={'size': 11})
	fig.savefig('TLinfty.png', bbox_inches="tight")
	plt.show()
	
	print("T( L=infty ) = ", TcInfty) 
	
	"""
	>python analyseDataP4.py Tinfty
	T( L=infty ) =  2.4244999999999997
	"""
	


