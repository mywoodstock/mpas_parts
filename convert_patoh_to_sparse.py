##
 #  Project:
 #
 #  File: convert_patoh_to_sparse.py
 #  Created: June 24, 2014
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 #
 #	Description: Convert PATOH hypergraph representation into sparse graph representation (row, col, val)
 ##

import sys, getopt

def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, "i:o:")
		if len(opts) < 2:
			raise getopt.GetoptError("Give arguments")
	except getopt.GetoptError:
		print "convert_patoh_to_sparse.py -i <inputfile> -o <outputfile>"
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-i"):
			inputfile = arg
		elif opt in ("-o"):
			outputfile = arg
		else:
			print "wrong argument"
			sys.exit(2)
	return inputfile, outputfile

# assuming the nets in input patoh file are ordered according to vertices (netid == vertexid)
def read_patoh_file(filename, data):
	ff = open(filename)
	head = ff.readline()
	header = ''
	while head.startswith("%"):
		header += head
		head = ff.readline()
	headwords = head.strip().split()
	if len(headwords) == 4:
		index = int(headwords[0])
		nnodes = int(headwords[1])
		nnets = int(headwords[2])
		nedges = int(headwords[3])
		if nnodes != nnets:
			print "error: need symmetric hypergraph (nnodes == nnets)"
			ff.close()
			sys.exit(2)
	else:
		print "wrong header in input patoh file"
		ff.close()
		sys.exit(2)
	count = 0
	while True:
		line = ff.readline()
		if not line: break
		count += 1
		net = map(int, line.strip().split())
		if index == 0:
			net = map(lambda x: x + 1, net)
		elif index == 1:
			net = net
		else:
			print 'error: wrong starting index specified in the header'
			ff.close()
			sys.exit(2)
		data[count] = net
	ff.close()
	if count != nnets:
		print 'error: mismatch in number of lines and nnets'
		sys.exit(2)
	return nnodes, nedges

## assuming symmetric data is already present
def convert_hyper_to_sparse(hyper_data, sparse_data):
	for net, nodes in hyper_data.iteritems():
		for node in nodes:
			sparse_data += [(net, node, 1)]


def write_sparse_file(filename, nnodes, nedges, sparse_data):
	ff = open(filename, 'w')
	head = '%%MatrixMarket matrix coordinate real general\n%\n'
	head += "%d\t%d\t%d\n" % (nnodes, nnodes, nedges)
	ff.write(head)
	count = 0
	for edge in sparse_data:
		record = '%d\t%d\t%f\n' % edge
		ff.write(record)
		count += 1
	ff.close()
	return count


infile, outfile = parse_arguments(sys.argv[1:])
hyper_data = {}
nnodes, nedges = read_patoh_file(infile, hyper_data)
sparse_data = []
npins = convert_hyper_to_sparse(hyper_data, sparse_data)
print "nodes = %d, edges = %d" % (nnodes, nedges)
write_sparse_file(outfile, nnodes, nedges, sparse_data)
