import matplotlib.pyplot as plt
import numpy as np

def relative_error(infile, color, label_, save_file = False, outfile = ''):
	"""
	Takes input file with log(relative error) for different # of steps
	Makes logplot over # steps and writes to outfile
	"""

	data = open(infile, "r").readlines()	
	error = np.zeros(len(data))
	
	for i in range(len(data)):
		error[i] = data[i]
	
	exponent = np.linspace(1,len(error),len(error))
	
	plt.plot(exponent, error, color, label=label_)
	if save_file:
		plt.xlabel('log(N)')
		plt.ylabel('log(relative error)')
		plt.legend()
		plt.savefig(outfile)

if __name__ == "__main__":
	
	relative_error("relativeError_slow.dat", 'g--', 'General algorithm')
	relative_error("relativeError_fast.dat", 'r:', 'Special algorithm')
	relative_error("relativeError_arma.dat", 'b--', 'Armadillo algorithm', True, "relative_error.png",)

	print("Your plot(s) were saved as:\n    relative_errors.png")
