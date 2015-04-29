import sys, getopt
import numpy as np

class Graph:

	graph = {}
	#A = 0.00002494
	#B = 0.0004956
	A = 2.494		## computation cost factor
	B = 49.56		## communication cost factor

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
			parts = map(int, line.strip().split())
			for part in parts:
				node += 1
				if part not in self.graph['partitions']:
					self.graph['partitions'][part] = {}
					self.graph['partitions'][part]['cells'] = set()
					self.graph['partitions'][part]['halo_cells'] = set()
				self.graph['partitions'][part]['cells'].add(node)
				self.graph['cells'][node]['partition'] = part
		pp.close()

	def get_halo_candidates(self, candidates, cell, layers):
		for neighbor in self.cells_on_cell(cell):
			candidates.add(neighbor)
			if layers > 0: self.get_halo_candidates(candidates, neighbor, layers - 1)

	def compute_halos(self, layers):
		print "computing halos ..."
		for part in range(0, self.num_partitions()):
			candidates = set()
			for cell in self.cells_in_partition(part):
				self.get_halo_candidates(candidates, cell, layers)
			self.graph['partitions'][part]['halo_cells'] = list(candidates - 
					set(self.graph['partitions'][part]['cells']))

	def print_cost(self, num_cells, num_halo_cells, label):
		#ccomp = (num_cells + num_halo_cells) * self.A * self.num_partitions()
		#ccomm = num_halo_cells * self.B * self.num_partitions()
		ccomp = (num_cells + num_halo_cells) * self.A
		ccomm = num_halo_cells * self.B
		print "**   computation cost " + label + ": " + str(ccomp)
		print "** communication cost " + label + ": " + str(ccomm)
		print "**         total cost " + label + ": " + str(ccomp + ccomm)

	def print_statistics(self):
		num_parts = self.num_partitions()
		temp_part_num_cells = []
		temp_part_num_halo_cells = []
		for part, partdata in self.graph['partitions'].iteritems():
			temp_part_num_cells.append(len(partdata['cells']))
			temp_part_num_halo_cells.append(len(partdata['halo_cells']))
		part_num_cells = np.array(temp_part_num_cells)
		part_num_halo_cells = np.array(temp_part_num_halo_cells)
		COMP = 1.0
		COMM = 1.0
		#part_total_num_cells = COMP * (part_num_cells + part_num_halo_cells) + COMM * part_num_halo_cells
		part_total_num_cells = part_num_cells + part_num_halo_cells
		print ("**            num cells:"
				+ " min = " + str(part_num_cells.min())
				+ " max = " + str(part_num_cells.max())
				+ " mean = " + str(part_num_cells.mean())
				+ " std-dev = "	+ str(part_num_cells.std())
				+ " imbalance = " + str(1.0 - float(part_num_cells.min())/part_num_cells.max()))
		print ("**       num halo cells:"
				+ " min = " + str(part_num_halo_cells.min())
				+ " max = " + str(part_num_halo_cells.max())
				+ " mean = " + str(part_num_halo_cells.mean())
				+ " std-dev = " + str(part_num_halo_cells.std())
				+ " imbalance = " + str(1.0 - float(part_num_halo_cells.min())/part_num_halo_cells.max()))
		print ("**      total num cells:"
				+ " min = " + str(part_total_num_cells.min())
				+ " max = " + str(part_total_num_cells.max())
				+ " mean = " + str(part_total_num_cells.mean())
				+ " std-dev = " + str(part_total_num_cells.std())
				+ " imbalance = " + str(1.0 - float(part_total_num_cells.min())/part_total_num_cells.max()))
		self.print_cost(part_num_cells.mean(), part_num_halo_cells.mean(), "mean")
		self.print_cost(part_num_cells.max(), part_num_halo_cells.max(), "max")
		self.print_cost(part_num_cells.min(), part_num_halo_cells.min(), "min")

	def print_detailed_statistics(self):
		print "** partitions data:"
		for part, pdata in self.graph['partitions'].iteritems():
			ncells = len(pdata['cells'])
			nhcells = len(pdata['halo_cells'])
			print ("    partition " + str(part) + " :: cells = " + str(ncells)
					+ " :: halo cells = " + str(nhcells)
					+ " :: total cells = " + str(nhcells + ncells))
#			print "    ==== cells: " + str(sorted(pdata['cells']))
#			print "    oooo halo cells: " + str(sorted(pdata['halo_cells']))
		self.print_statistics()

	def process(self):
		self.compute_halos(3)
		self.print_detailed_statistics()
	
	def printall(self):
		print self.graph


def parse_arguments(argv):
	try:
		opts, args, = getopt.getopt(argv, "g:p:")
		if len(opts) < 2:
			raise getopt.GetoptError("Not enough arguments provided")
	except getopt.GetoptError:
		print "construct_graph_data.py -g <graph_file> -p <partition_file>"
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

