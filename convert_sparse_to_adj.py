##
 #  Project:
 #
 #  File: convert_sparse_to_adj.py
 #  Created: Mar 08, 2014
 #  Modified: Mon 17 Mar 2014 09:56:29 AM PDT
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

import sys, getopt

def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, "i:o:")
		if len(opts) < 2:
			raise getopt.GetoptError("Give arguments")
	except getopt.GetoptError:
		print "convert_sparse_to_adj.py -i <inputfile> -o <outputfile>"
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-i"):
			inputfile = arg
		elif opt in ("-o"):
			outputfile = arg
	return inputfile, outputfile


def read_sparse_file(filename, data):
	ff = open(filename)
	head = ff.readline()		# this has num nodes and num edges
	header = ""
	while head.startswith("%"):
		header += head
		head = ff.readline()
	headwords = head.strip().split()
	if len(headwords) == 2:
		nodes = int(headwords[0])
		edges = int(headwords[1])
	elif len(headwords) == 3:
		nodes = int(headwords[0])
		nodes2 = int(headwords[1])
		if nodes != nodes2:
			print "error: input matrix is not square"
			sys.exit(2)
		edges = int(headwords[2])
	count = 0
	while True:
		line = ff.readline()
		if not line: break
		count += 1
		words = line.strip().split()
		data += [(int(words[0]), int(words[1]), int(words[2]))]
	if count != edges:
		print "error: mismatching number of edges and data in the file: %d, %d" % (edges, count)
		ff.close()
		sys.exit(2)
	ff.close()
	return nodes, edges


def convert_sparse_to_adj(sparse_data, adj_data):
	for row, col, val in sparse_data:
		if row < col:
			if row not in adj_data:
				adj_data[row] = []
			adj_data[row] += [col]
			# symmetrize
			if col not in adj_data:
				adj_data[col] = []
			adj_data[col] += [row]


def write_adj_file(filename, nodes, edges, adj_data):
	ff = open(filename, 'w')
	head = "%d\t%d\n" % (nodes, edges)
	ff.write(head)
	count = 1
	for row, neighbors in adj_data.items():
		while row > count:
			record = "\n"
			ff.write(record)
			count += 1
		record = ""
		for col in neighbors:
			record += "%d\t" % col
		record += "\n"
		ff.write(record)
		count += 1
	ff.close()


infile, outfile = parse_arguments(sys.argv[1:])
sparse_data = []
nnodes, nedges = read_sparse_file(infile, sparse_data)
print "nodes = %d, edges = %d" % (nnodes, nedges)
adj_data = {}
convert_sparse_to_adj(sparse_data, adj_data)
write_adj_file(outfile, nnodes, nedges, adj_data)
