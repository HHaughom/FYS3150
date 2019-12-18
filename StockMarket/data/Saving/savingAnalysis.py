# Python script for data analysing

import math as math
import numpy as np
import matplotlib.pyplot as plt

def w(m, m0):
    beta = 1.0/m0;
    return beta*np.exp(-beta*m)

def n(saving):
    return 1 + 3*saving/(1-saving)

def a(n):
    return n**n/math.gamma(n)

def Pn(x, saving):
    var = n(saving)
    return a(var)*x**(var-1)*np.exp(-var*x)


if __name__ == '__main__':

    lambda0 = np.loadtxt("lambda0.dat")
    lambda025 = np.loadtxt("lambda025.dat")
    lambda05 = np.loadtxt("lambda05.dat")
    lambda09 = np.loadtxt("lambda09.dat")
    x = np.linspace(0,10,len(lambda0))

    plt.figure(1)
    plt.scatter(x, lambda0, marker=".", label="lambda = 0")
    plt.scatter(x, lambda025, marker="+", label="lambda = 0.25")
    plt.scatter(x, lambda05, marker="x", label="lambda = 0.5")
    plt.scatter(x, lambda09, marker="1", label="lambda = 0.9")


    y1 = w(x,1.0)
    y2 = Pn(x, 0.25)
    y3 = Pn(x, 0.5)
    y4 = Pn(x, 0.9)

    sum1 = 0
    sum2 = 0
    sum3 = 0
    sum4 = 0
    for i in range(len(x)):
        sum1 += y1[i]
        sum2 += y2[i]
        sum3 += y3[i]
        sum4 += y4[i]

    y1 /= sum1
    y2 /= sum2
    y3 /= sum3
    y4 /= sum4

    plt.plot(x,y1)
    plt.plot(x,y2)
    plt.plot(x,y3)
    plt.plot(x,y4)
    plt.legend(loc="best")
    plt.xlabel('m')
    plt.ylabel('f(m)')
    plt.show()

    # Fka greie
    plt.figure(2)
    plt.scatter(x, lambda0, label="lambda0")
    plt.yscale("log")
    plt.xlabel('m')
    plt.ylabel('log(f(m))')
    plt.show()
