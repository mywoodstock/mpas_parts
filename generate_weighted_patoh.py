##
 #  Project:
 #
 #  File: generate_weighted_patoh.py
 #  Created: Mar 08, 2014
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
		print "convert_patoh_to_mpas.py -i <inputfile> -o <outputfile>"
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-i"):
			inputfile = arg
		elif opt in ("-o"):
			outputfile = arg
	return inputfile, outputfile


def read_hypergraph(filename, data):
	ff = open(filename)
	line = ff.readline()		## first line has [ base, #nodes, #nets, #pins, weightscheme=1 ]
	if not line: return
	words = map(int, line.strip().split())
	baseidx = words[0]
	nnodes = words[1]
	nnets = words[2]
	npins = words[3]
	if len(words) == 5: wscheme = words[4]
	else: wscheme = 0
	count = 0
	while True:
		line = ff.readline()
		if not line: break
		if count == nnets: break	## ignore the last weights line
		pins = map(lambda x: x - baseidx, map(int, line.strip().split()))
		data[count] = pins	 
		count += 1
	ff.close()
	if count != nnets:
		print 'error: mismatch in number of nets and lines'
		sys.exit(2)
	return nnodes, nnets, npins


def write_weighted_hypergraph(filename, data, nnodes, nnets, npins):
	if nnodes != nnets: print 'warning: nnodes and nnets do not match: %d, %d\n' % (nnodes, nnets)
	weights = []
	ff = open(filename, 'w')
	header = '0 %d %d %d 1\n' % (nnodes, nnets, npins)
	ff.write(header)
	for net, pins in data.iteritems():
		weights.append(len(pins))
		record = ''
		for p in pins: record += str(p) + ' '
		record += '\n'
		ff.write(record)
	record = ''
	for w in weights: record += str(w) + ' '
	record += '\n'
	ff.write(record)
	ff.close()


infile, outfile = parse_arguments(sys.argv[1:])
data = {}	## mapping { net -> list of nodes in net }
## the obtained data always has base index as 0
nnodes, nnets, npins = read_hypergraph(infile, data)
write_weighted_hypergraph(outfile, data, nnodes, nnets, npins)
