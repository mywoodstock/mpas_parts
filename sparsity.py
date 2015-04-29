from matplotlib import pyplot as plt
from scipy.io import mmread
import os.path
import sys

def showSparsity(matFile):
	M = mmread(matFile)
	plt.title(matFile)
	plt.spy(M, precision=1e-3, marker='.', markersize=5)
	plt.show()

if __name__ == "__main__":
	showSparsity(sys.argv[1])
