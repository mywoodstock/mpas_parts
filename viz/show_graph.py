##
 #  Project:
 #
 #  File: show_graph.py
 #  Created: Mar 09, 2014
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

import os
import networkx as nx
print "Display is: %s" % (os.environ['DISPLAY'],)
import matplotlib
print "Default backend is: %s" % (matplotlib.get_backend(),)
#matplotlib.use("gtk")
#print "Backend is now: %s" % (matplotlib.get_backend(),)

from scipy import io
from matplotlib import pyplot, patches

import sys, getopt

def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, "hi:")
		if len(opts) < 1:
			raise getopt.GetoptError("Give arguments")
	except getopt.GetoptError:
		print 'show_graph.py -i <inputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'show_graph.py -i <inputfile>'
			sys.exit()
		elif opt in ("-i"):
			inputfile = arg
	return inputfile


def draw_adjacency_matrix(G, node_order = None, partitions = [], colors = []):
	"""
	- G is a netorkx graph
	- node_order (optional) is a list of nodes, where each node in G appears exactly once
	- partitions is a list of node lists, where each node in G appears in exactly one node list
	- colors is a list of strings indicating what color each partition should be
	If partitions is specified, the same number of colors needs to be specified.
	"""
	adjacency_matrix = nx.to_numpy_matrix(G, dtype=bool, nodelist=node_order)
	#Plot adjacency matrix in toned-down black and white
	fig = pyplot.figure(figsize = (5, 5)) # in inches
	pyplot.imshow(adjacency_matrix, cmap = "Greys", interpolation = "none")
	# The rest is just if you have sorted nodes by a partition and want to
	# highlight the module boundaries
	assert len(partitions) == len(colors)
	ax = pyplot.gca()
	for partition, color in zip(partitions, colors):
		current_idx = 0
		for module in partition:
			ax.add_patch(patches.Rectangle((current_idx, current_idx),
										  len(module), # Width
										  len(module), # Height
										  facecolor="none",
										  edgecolor=color,
										  linewidth="1"))
			current_idx += len(module)


infile = parse_arguments(sys.argv[1:])
A = io.mmread(infile)
G = nx.from_scipy_sparse_matrix(A)
draw_adjacency_matrix(G)
pyplot.show()
