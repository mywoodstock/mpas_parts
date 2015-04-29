import sys, getopt
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

def parse_arguments(argv):
	try:
		opts, args, = getopt.getopt(argv, "i:o:")
		if len(opts) < 2:
			raise getopt.GetoptError("Not enough arguments provided")
	except getopt.GetoptError:
		print "count_part_cells.py -i <part_file> -o <output_plot_file>"
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-i"):
			inputfile = arg
		if opt in ("-o"):
			outfile = arg
	return inputfile, outfile


def read_parts(f, parts):
	ff = open(f)
	while True:
		line = ff.readline()
		if not line: break;
		parts += [ int(line.strip().split()[0]) ]
	ff.close()


def count_cells(f, counts):
	ff = open(f)
	while True:
		line = ff.readline()
		if not line: break
		part = int(line.strip().split()[0])
		if part not in counts: counts[part] = 0
		counts[part] += 1
	ff.close()


def plot_part_histogram(outfile, data):
	data = np.array(data)
	counts,_,_ = plt.hist(data, bins=(max(data) - min(data) + 1), rwidth=0.75, color='g', alpha=0.75, lw=0)
	mean = counts.mean()
	std = counts.std()
	print mean, std
	x = range(0, max(data) - min(data) + 1)
	y = np.zeros(len(x))
	y = y + mean
	plt.plot(x, y, lw = 1, alpha = 0.5, color = 'k')
	plt.xlabel("Partition Number")
	plt.ylabel("nCells")
	plt.xlim(0, max(data) - min(data))
	plt.ylim(0, plt.ylim()[1] + 50)
	plt.gcf().set_size_inches(15, 6)
	plt.savefig(outfile)


infile, outfile = parse_arguments(sys.argv[1:])
#part_counts = {}
parts = []
#count_cells(infile, part_counts)
read_parts(infile, parts)
plot_part_histogram(outfile, parts)


