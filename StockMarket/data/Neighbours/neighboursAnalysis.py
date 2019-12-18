# Analysis script for Neighbours data
import numpy as np
import matplotlib.pyplot as plt

#N = 1000
#lambda = 0
#alpha = 0, 0.5, 1.0, 1.5, 2.0

N1000alfa0lambda0 = np.loadtxt("N1000alfa0lambda0.dat")
N1000alfa05lambda0 = np.loadtxt("N1000alfa05lambda0.dat")
N1000alfa1lambda0 = np.loadtxt("N1000alfa1lambda0.dat")
N1000alfa15lambda0 = np.loadtxt("N1000alfa15lambda0.dat")
N1000alfa2lambda0 = np.loadtxt("N1000alfa2lambda0.dat")

N1000lambda0 = [N1000alfa05lambda0, N1000alfa1lambda0, N1000alfa15lambda0, N1000alfa2lambda0]

#N = 1000
#lambda = 0.5
#alpha = 0.5, 1.0, 1.5, 2.0

N1000alfa05lambda05 = np.loadtxt("N1000alfa05lambda05.dat")
N1000alfa1lambda05 = np.loadtxt("N1000alfa1lambda05.dat")
N1000alfa15lambda05 = np.loadtxt("N1000alfa15lambda05.dat")
N1000alfa2lambda05 = np.loadtxt("N1000alfa2lambda05.dat")

N1000lambda05 = [N1000alfa05lambda05, N1000alfa1lambda0, N1000alfa15lambda0, N1000alfa2lambda0]

#N = 500
#lambda = 0
#alpha = 0.5, 1.0, 1.5, 2.0

N500alfa05lambda0 = np.loadtxt("N500alfa05lambda0.dat")
N500alfa1lambda0 = np.loadtxt("N500alfa1lambda0.dat")
N500alfa15lambda0 = np.loadtxt("N500alfa15lambda0.dat")
N500alfa2lambda0 = np.loadtxt("N500alfa2lambda0.dat")

N500lambda0 = [N500alfa05lambda0, N500alfa1lambda0, N500alfa15lambda0, N500alfa2lambda0]
x = np.linspace(0, 20, len(N1000alfa0lambda0))
labels = ["alpha = 0.5", "alpha = 1.0", "alpha = 1.5", "alpha = 2.0"]


if __name__ == '__main__':
    for para in [N1000lambda0, N1000lambda05, N500lambda0]:
    	for data in range(len(para)):
    		plt.plot(x, para[data], ':', label = labels[data])

    	plt.yscale("log")
    	plt.xscale("log")
    	plt.legend(loc="best")
    	plt.xlabel('m')
    	plt.ylabel('f(m)')
    	plt.show()
