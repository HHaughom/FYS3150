import matplotlib.pyplot as plt
import numpy as np
import sys

def u_plot(filename, outfile, flipped=False):
	"""
	Function takes a datafile with solution u and a figure filename as input.
	Plots numeric and analytic solution u and u_analytical.
	Writes plot to figure file. 
	"""

	data = open(filename, "r").readlines()
	u = np.zeros(len(data))
	if flipped == True:
		for i in range(len(u)):
			u[i] = data[-i-1]
	else:
		for i in range(len(u)):
			u[i] = data[i]
	x = np.linspace(0,1,len(u))
	u_analytical = 1-(1-np.exp(-10))*x-np.exp(-10*x)
	
	plt.plot(x,u_analytical, 'g-', label='Analytic')
	plt.plot(x,u, 'r:', label='Numeric')
	plt.legend()
	plt.xlabel('x'); plt.ylabel('u')
	plt.savefig(outfile)
	plt.clf()
	

if __name__ == "__main__":

	#Initialize error to store plots 
	plots_created = []
	
	#Asking user to choose algorithm 
	succsess = False
	while succsess == False:
		print("What algorithm do you want?\n1) General algorithm\n2) Specialized algorithm\n3) Armadillo LU decomposition")
		algo = int(input("I want: "))

		if algo in range(1,4):
			print("\nAlrighty. Good choice.")
			succsess = True
		else:
			print("You need to choose an algorithm between 1 and 3.")	

	succsess = False

	#Asking user to choose number of integration steps
	#Chice between 1e1 and 1e7 or all of them 
	while succsess == False:
		print("\nChoose number of integration steps.\nPick an order of magnitude between 1 and 7.\nOr type \"all\" to plot all magnitudes.")
		exponent = input("Order of magnitude: ")
	
		if (exponent == "all" or int(exponent) in range(1,8)) :
			succsess = True
		else:
			print("\nYou need to choose number of integration steps. Try again.")

	if algo == 1:

		if exponent == "all":
	
			for exp in range(1,8):
				n = int(10**exp)
				u_plot("u_data_slow" + str(n) + "n.dat", "slow_solution" + str(n)+ "n.png", True)
				plots_created.append(("slow_solution" + str(n)+ "n.png"))
		else:
			n = int(10**int(exponent))
			u_plot("u_data_slow" + str(n) + "n.dat", "slow_solution" + str(n)+ "n.png", True)
			plots_created.append(("slow_solution" + str(n)+ "n.png"))
			
	elif algo == 2:

		if exponent == "all":

			for exp in range(1,8):
				n = int(10**exp)
				u_plot("u_data_fast" + str(n) + "n.dat", "fast_solution" + str(n)+ "n.png", True)
				plots_created.append(("fast_solution" + str(n)+ "n.png"))

		else:
			n = int(10**int(exponent))
			u_plot("u_data_fast" + str(n) + "n.dat", "fast_solution" + str(n)+ "n.png", True)
			plots_created.append(("fast_solution" + str(n)+ "n.png"))
		
	elif algo == 3:

		if exponent == "all":

			for exp in range(1,5):
				n = int(10**exp)
				u_plot("u_data_arma" + str(n) + "n.dat", "arma_solution" + str(n)+ "n.png")
				plots_created.append(("arma_solution" + str(n)+ "n.png"))

		else:
			n = int(10**int(exponent))
			u_plot("u_data_arma" + str(n) + "n.dat", "arma_solution" + str(n)+ "n.png")
			plots_created.append(("arma_solution" + str(n)+ "n.png"))

	print("Your plot(s) were saved as:\n")
	for plot in plots_created:
		print("    " + plot)
