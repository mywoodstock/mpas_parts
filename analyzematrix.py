##
 #  Project:
 #
 #  File: analyzematrix.py
 #  Created: Apr 28, 2014
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

import sys, getopt
import numpy as np
from sets import Set

class Graph:

	graph = {}

	def __init__(self):
		self.graph['partitions'] = {}
		self.graph['cells'] = {}

	def num_cells(self):
		return len(self.graph['cells'])

	def cells_on_cell(self, cell):
		if cell not in self.graph['cells']: return []
		return self.graph['cells'][cell]['neighbors']
	
	def cell_partition(self, cell):
		if cell not in self.graph['cells']: return -1
		return self.graph['cells'][cell]['partition']

	def cells_in_partition(self, part):
		if part not in self.graph['partitions']: return []
		return self.graph['partitions'][part]['cells']
	
	def num_partitions(self):
		return len(self.graph['partitions'])

	def num_cells_in_partition(self, part):
		if part not in self.graph['partitions']: return 0
		return len(self.graph['partitions'][part]['cells'])

	def read_graph(self, gfile):
		print "reading graph ... "
		gg = open(gfile)
		line = gg.readline()	# first line has num cells and num edges
		words = line.strip().split()
		num_cells = int(words[0])
		num_edges = int(words[1])
		node = 0
		while True:
			line = gg.readline()
			if not line: break
			node += 1
			words = line.strip().split()
			if len(words) > 0: neighbors = map(int, words)
			else: neighbors = []
			self.graph['cells'][node] = {}
			self.graph['cells'][node]['neighbors'] = neighbors
			self.graph['cells'][node]['partition'] = 0			# initialize
			for neighbor in neighbors:
				if neighbor not in self.graph['cells']: self.graph['cells'][neighbor] = {}
				if 'neighbors' not in self.graph['cells'][neighbor]:
					self.graph['cells'][neighbor]['neighbors'] = []
				if node not in self.graph['cells'][neighbor]['neighbors']:
					self.graph['cells'][neighbor]['neighbors'] += [ node ]
		gg.close()
		if node != num_cells:
			print "number of cells mismatch in graph file"
			sys.exit(2)

	def read_partitions(self, pfile):
		print "reading partitions ... "
		pp = open(pfile)
		node = 0
		while True:
			line = pp.readline()
			if not line: break
			node += 1
			part = int(line.strip().split()[0])
			if part not in self.graph['partitions']:
				self.graph['partitions'][part] = {}
				self.graph['partitions'][part]['cells'] = set()
			self.graph['partitions'][part]['cells'].add(node)
			self.graph['cells'][node]['partition'] = part
		pp.close()


	def calculate_partition_off_diagonals(self, part, pdata):
		neighbors = Set()
		for cell in pdata['cells']:
			neighbors.update(self.graph['cells'][cell]['neighbors'])
		return list(neighbors - Set(pdata['cells']))


	def calculate_off_diagonals(self):
		print "calculating number of off-diagonal nodes ..."
		## find min and max node ids for each partition
		for part, pdata in self.graph['partitions'].iteritems():
			pmin = min(pdata['cells'])
			pmax = max(pdata['cells'])
			pdata['offd_cells'] = self.calculate_partition_off_diagonals(part, pdata)


	def print_statistics(self):
		num_parts = self.num_partitions()
		temp_part_num_cells = []
		temp_part_num_offd_cells = []
		for part, partdata in self.graph['partitions'].iteritems():
			temp_part_num_cells.append(len(partdata['cells']))
			temp_part_num_offd_cells.append(len(partdata['offd_cells']))
		part_num_cells = np.array(temp_part_num_cells)
		part_num_offd_cells = np.array(temp_part_num_offd_cells)
		part_total_num_cells = part_num_cells + part_num_offd_cells
		print ("**            num cells:"
				+ " min = " + str(part_num_cells.min())
				+ " max = " + str(part_num_cells.max())
				+ " mean = " + str(part_num_cells.mean())
				+ " std-dev = "	+ str(part_num_cells.std())
				+ " imbalance = " + str(1.0 - float(part_num_cells.min())/part_num_cells.max()))
		print ("**       num offd cells:"
				+ " min = " + str(part_num_offd_cells.min())
				+ " max = " + str(part_num_offd_cells.max())
				+ " mean = " + str(part_num_offd_cells.mean())
				+ " std-dev = " + str(part_num_offd_cells.std())
				+ " imbalance = " + str(1.0 - float(part_num_offd_cells.min())/part_num_offd_cells.max()))
		print ("**      total num cells:"
				+ " min = " + str(part_total_num_cells.min())
				+ " max = " + str(part_total_num_cells.max())
				+ " mean = " + str(part_total_num_cells.mean())
				+ " std-dev = " + str(part_total_num_cells.std())
				+ " imbalance = " + str(1.0 - float(part_total_num_cells.min())/part_total_num_cells.max()))

	def print_detailed_statistics(self):
		print "** partitions data:"
		for part, pdata in self.graph['partitions'].iteritems():
			ncells = len(pdata['cells'])
			nocells = len(pdata['offd_cells'])
			print ("    partition " + str(part) + " :: cells = " + str(ncells)
					+ " :: offd cells = " + str(nocells)
					+ " :: total cells = " + str(nocells + ncells)
					+ " :: comm/comp ratio = " + str(float(nocells) / (ncells + nocells)))
#			print "    ==== cells: " + str(pdata['cells'])
#			print "    oooo offd cells: " + str(pdata['offd_cells'])
		self.print_statistics()

	def process(self):
		self.calculate_off_diagonals()
		self.print_detailed_statistics()
	
	def printall(self):
		print self.graph


def parse_arguments(argv):
	try:
		opts, args, = getopt.getopt(argv, "g:p:")
		if len(opts) < 2:
			raise getopt.GetoptError("Not enough arguments provided")
	except getopt.GetoptError:
		print "analyzematrix.py -g <graph_file> -p <partition_file>"
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-g"):
			graphfile = arg
		if opt in ("-p"):
			partfile = arg
	return graphfile, partfile



if __name__ == "__main__":
	gfile, pfile = parse_arguments(sys.argv[1:])
	graph = Graph()
	graph.read_graph(gfile)
	graph.read_partitions(pfile)
	graph.process()

