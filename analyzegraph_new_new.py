import sys, getopt
import numpy as np

import networkx as nx
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


## colors
myblue		= "#3355A9"		# blue
myred		= "#bc2318"		# red
mygreen		= "#4d9b59"		# green
myyellow	= "#d2a41a"		# yellow
mylblue		= "#7799ff"		# blue
mylred		= "#ff675c"		# red
mylgreen	= "#8fdf9d"		# green
mygray		= "#999999"		# gray

## bruegel
bruegel01	= "#d6553d"
bruegel02	= "#755e24"
bruegel03	= "#514639"
bruegel04	= "#986b5c"
bruegel05	= "#a5c596"
bruegel06	= "#749867"
bruegel07	= "#bbae76"
bruegel08	= "#fff2c2"
bruegel09	= "#e8df96"
bruegel10	= "#909071"

## raphael
raphael01	= "#76a79f"
raphael02	= "#e6bfa9"
raphael03	= "#f3e6bd"
raphael04	= "#f4e7cf"
raphael05	= "#8e3649"
raphael06	= "#7e8ba9"
raphael07	= "#70a56c"
raphael08	= "#b76669"
raphael09	= "#546288"
raphael10	= "#446224"
raphael11	= "#bfc0c1"

## motherwell
motherwell01	= "#928179"
motherwell02	= "#6a9cad"
motherwell03	= "#c19f7f"
motherwell04	= "#b1033d"
motherwell05	= "#8f9061"
motherwell06	= "#a5a677"
motherwell07	= "#7e6b71"
motherwell08	= "#7b7b6b"

## colormaps
motherwelldict1 = {
			'red':		((0.0, 0.573, 0.573),
						(1.0, 0.416, 0.416)),
			'green':	((0.0, 0.506, 0.506),
						(1.0, 0.612, 0.612)),
			'blue':		((0.0, 0.475, 0.475),
						(1.0, 0.678, 0.678))
			}
motherwellcmap1 = LinearSegmentedColormap('motherwellcmap1', motherwelldict1)

motherwelldict2 = {
			'red':		((0.0, 0.573, 0.573),
						(1.0, 0.694, 0.694)),
			'green':	((0.0, 0.506, 0.506),
						(1.0, 0.012, 0.012)),
			'blue':		((0.0, 0.475, 0.475),
						(1.0, 0.239, 0.239))
			}
motherwellcmap2 = LinearSegmentedColormap('motherwellcmap2', motherwelldict2)

motherwelldict3 = {
			'red':		((0.0, 0.416, 0.416),
						(1.0, 0.694, 0.694)),
			'green':	((0.0, 0.612, 0.612),
						(1.0, 0.012, 0.012)),
			'blue':		((0.0, 0.678, 0.678),
						(1.0, 0.239, 0.239))
			}
motherwellcmap3 = LinearSegmentedColormap('motherwellcmap3', motherwelldict3)

motherwelldict4 = {
			'red':		((0.0, 0.757, 0.757),
						(1.0, 0.694, 0.694)),
			'green':	((0.0, 0.624, 0.624),
						(1.0, 0.012, 0.012)),
			'blue':		((0.0, 0.498, 0.498),
						(1.0, 0.239, 0.239))
			}
motherwellcmap4 = LinearSegmentedColormap('motherwellcmap4', motherwelldict4)


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
			parts = map(int, line.strip().split())
			for part in parts:
				node += 1
				if part not in self.graph['partitions']:
					self.graph['partitions'][part] = {}
					self.graph['partitions'][part]['cells'] = set()
					self.graph['partitions'][part]['halo_cells'] = set()
					self.graph['partitions'][part]['neighbors'] = set()
				self.graph['partitions'][part]['cells'].add(node)
				self.graph['cells'][node]['partition'] = part
		pp.close()

	def get_halo_candidates(self, candidates, cell, layers):
		for neighbor in self.cells_on_cell(cell):
			candidates.add(neighbor)
			if layers > 0: self.get_halo_candidates(candidates, neighbor, layers - 1)

	## compute number of halo cells using recursion for 'layers' layers of halos
	def compute_halos(self, layers):
		print "computing halos ..."
		for part in range(0, self.num_partitions()):
			candidates = set()
			for cell in self.cells_in_partition(part):
				self.get_halo_candidates(candidates, cell, layers)
			self.graph['partitions'][part]['halo_cells'] = list(candidates - 
					set(self.graph['partitions'][part]['cells']))
	
	def compute_part_neighbors(self):
		nodes = self.graph['cells']
		for part, pdata in self.graph['partitions'].iteritems():
			for hcell in pdata['halo_cells']:
				pdata['neighbors'].add(nodes[hcell]['partition'])

	def print_statistics(self):
		num_parts = self.num_partitions()
		temp_part_num_cells = []
		temp_part_num_halo_cells = []
		temp_part_num_neighbors = []
		for part, partdata in self.graph['partitions'].iteritems():
			temp_part_num_cells.append(len(partdata['cells']))
			temp_part_num_halo_cells.append(len(partdata['halo_cells']))
			temp_part_num_neighbors.append(len(partdata['neighbors']))
		part_num_cells = np.array(temp_part_num_cells)
		part_num_halo_cells = np.array(temp_part_num_halo_cells)
		part_num_neighbors = np.array(temp_part_num_neighbors)
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
		print ("**            neighbors:"
				+ " min = " + str(part_num_neighbors.min())
				+ " max = " + str(part_num_neighbors.max())
				+ " mean = " + str(part_num_neighbors.mean())
				+ " std-dev = " + str(part_num_neighbors.std())
				+ " imbalance = " + str(1.0 - float(part_num_neighbors.min())/part_num_neighbors.max()))
		nparts = len(self.graph['partitions'])

	def print_detailed_statistics(self):
		print "** partitions data:"
		for part, pdata in self.graph['partitions'].iteritems():
			ncells = len(pdata['cells'])
			nhcells = len(pdata['halo_cells'])
			nneighbors = len(pdata['neighbors'])
			print ("    partition " + str(part) + " :: cells = " + str(ncells)
					+ " :: halo cells = " + str(nhcells)
					+ " :: total cells = " + str(nhcells + ncells)
					+ " :: neighbors = " + str(nneighbors))
		self.print_statistics()

	def printall(self):
		print self.graph
	
	def process(self):
		self.compute_halos(3)
		self.compute_part_neighbors()
		self.print_detailed_statistics()
	
	def generate_partition_network_graph(self, oprefix):
		print 'generating visuals ...'
		pgraph = nx.Graph()
		pgraph.add_nodes_from(np.arange(0, len(self.graph['partitions'])))
		for part, pdata in self.graph['partitions'].iteritems():
			pgraph.add_edges_from(map(lambda x: (part,x), pdata['neighbors']))
		# x = int(np.sqrt(pgraph.number_of_nodes()))
		# y = pgraph.number_of_nodes() / x
		# h = nx.empty_graph(24)
		# pos = nx.fruchterman_reingold_layout(pgraph)
		pos = nx.graphviz_layout(pgraph, prog='neato') 		## prog = twopi, gvcolor, wc, ccomps, tred, sccmap, fdp, circo, neato, acyclic, nop, gvpr, dot, sfdp
		plt.clf()
		node_sizes = [ len(vdata['cells'])/5. for v, vdata in self.graph['partitions'].iteritems() ]
		node_colors = [ len(vdata['cells'])/5. for v, vdata in self.graph['partitions'].iteritems() ]
		line_widths = [ len(vdata['halo_cells'])/200. for v, vdata in self.graph['partitions'].iteritems() ]
		# nx.draw(pgraph, pos, cmap=motherwellcmap4, node_size=node_sizes, linewidths=line_widths, node_color=node_colors, width=1, edge_color=motherwell01, with_labels=False)
		edges = nx.draw_networkx_edges(pgraph, pos, width=1, edge_color=motherwell01, alpha=0.5, label=None)
		nodes = nx.draw_networkx_nodes(pgraph, pos, cmap=motherwellcmap4, node_size=node_sizes, linewidths=line_widths, node_color=node_colors, width=1, edge_color=motherwell01, with_labels=False)
		nodes.set_edgecolor('#444444')
		plt.gcf().set_size_inches(7,7)
		plt.axis('off')
		plt.savefig(oprefix+'_pgraph.pdf', bbox_inches='tight')
		
	def generate_grid_network_graph(self, oprefix):
		return
	
	def save_partition_info(self, oprefix):
		## already have graph.info and partitions
		## save halos and partition neighbor graph
		hfilename = oprefix + '_partition_halos.info'
		ff = open(hfilename, 'w')
		p = str(self.num_partitions()) + '\n'
		ff.write(p)
		for part, pdata in sorted(self.graph['partitions'].iteritems()):
			phalo = ' '.join(map(str, pdata['halo_cells'])) + '\n'
			ff.write(phalo)
		ff.close()
		pfilename = oprefix + '_partition_graph.info'
		ff = open(pfilename, 'w')
		p = str(self.num_partitions()) + '\n'
		ff.write(p)
		for part, pdata in sorted(self.graph['partitions'].iteritems()):
			pneighbors = ' '.join(map(str, pdata['neighbors'])) + '\n'
			ff.write(pneighbors)
		ff.close()
	
	def read_partition_info(self, oprefix):
		print "reading halos ..."
		hfilename = oprefix + '_partition_halos.info'
		ff = open(hfilename, 'r')
		nparts = int(ff.readline().strip().split()[0])
		part = 0
		while True:
			line = ff.readline()
			if not line: break
			self.graph['partitions'][part]['halo_cells'] = map(int, line.strip().split())
			part += 1
		ff.close()
		if nparts != part:
			print 'error: mismatching number of partitions'
			sys.exit(2)
		print "reading partition graph ..."
		pfilename = oprefix + '_partition_graph.info'
		ff = open(pfilename, 'r')
		nparts = int(ff.readline().strip().split()[0])
		part = 0
		while True:
			line = ff.readline()
			if not line: break
			self.graph['partitions'][part]['neighbors'] = map(int, line.strip().split())
			part += 1
		ff.close()
		if nparts != part:
			print 'error: mismatching number of partitions'
			sys.exit(2)



def parse_arguments(argv):
	try:
		opts, args, = getopt.getopt(argv, "g:p:o:r:")
		if len(opts) < 4:
			raise getopt.GetoptError("Not enough arguments provided")
	except getopt.GetoptError:
		print "model_new.py -g <graph_file> -p <partition_file> -o <output_prefix> -r <read or compute>"
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-g"):
			graphfile = arg
		if opt in ("-p"):
			partfile = arg
		if opt in ("-o"):
			oprefix = arg
		if opt in ("-r"):
			read = arg
	return graphfile, partfile, oprefix, read



if __name__ == "__main__":
	gfile, pfile, oprefix, read = parse_arguments(sys.argv[1:])
	graph = Graph()
	graph.read_graph(gfile)
	graph.read_partitions(pfile)
	if read == '1':
		graph.read_partition_info(oprefix)
	else:
		graph.process()
		graph.save_partition_info(oprefix)
	graph.generate_partition_network_graph(oprefix)
	graph.generate_grid_network_graph(oprefix)
