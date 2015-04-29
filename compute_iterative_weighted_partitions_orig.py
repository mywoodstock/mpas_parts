##
 #  Project:
 #
 #  File: generate_weighted_patoh.py
 #  Created: Mar 08, 2014
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

import sys, getopt
import subprocess
import numpy as np
import pandas as pd
import random

from sets import Set
from operator import add, mul


A = 0.546	## computeation cost factor obtained from performance model
B = 0.78124	## communication cost factor obtained from performance model
TEMPERATURE = 1e-2


def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, 'i:p:o:n:t:')
		if len(opts) < 5 or len(opts) > 5:
			raise getopt.GetoptError('Give arguments')
	except getopt.GetoptError:
		print 'compute_iterative_weighted_partitions.py -i <hypergraphfile> -p <partitionfile> -o <outputprefix> -n <nparts> -t <niterations>'
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-i'): inputfile = arg
		elif opt in ('-o'): outputfile = arg
		elif opt in ('-p'): partfile = arg
		elif opt in ('-n'): nparts = int(arg)
		elif opt in ('-t'): niter = int(arg)
	return inputfile, partfile, outputfile, nparts, niter


def read_partitioned_hypergraph(filename, partfile, data, parts):
	## first read the hypergraph
	ff = open(filename)
	line = ff.readline()		## first line has [ base, #nodes, #nets, #pins, weightscheme=1 ]
	if not line: return
	words = map(int, line.strip().split())
	if len(words) < 4 or len(words) > 6:
		print 'error: invalid header in hypergraph file'
		sys.exit(2)
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

	## and then read the partition file
	count = 0
	ff = open(partfile)
	while True:
		line = ff.readline()
		if not line: break
		words = map(int, line.strip().split())
		for word in words:
			parts[count] = word
			count += 1
	ff.close()
	return nnodes, nnets, npins


def read_partitioning(partfile, parts):
	count = 0
	ff = open(partfile)
	while True:
		line = ff.readline()
		if not line: break
		words = map(int, line.strip().split())
		for word in words:
			parts[count] = word
			count += 1
	ff.close()


def compute_node_weights(data, parts, weights):
	# assuming one-to-one mapping nets <-> nodes
	for net, pins in data.iteritems():
		mypart = parts[net]
		for pin in pins:
			if parts[pin] != mypart: weights[pin] += 1
#	print "** new weights: ", weights
	return


def calculate_parts_weights(nparts, parts, weights, pweights):
	for node in range(0, len(weights)):
		part = parts[node]
		if part not in pweights: pweights[part] = 0
		pweights[part] += weights[node]
	for p in range(0, nparts):
		if p not in pweights: pweights[p] = 0


def construct_parts_halos(data, nparts, parts, halos):
	for net, pins in data.iteritems():		## pins represent all the nodes in the halos
		mypart = parts[net]
		if mypart not in halos: halos[mypart] = Set()
		for pin in pins:
			if parts[pin] != mypart: halos[mypart].add(pin)
	for p in range(0, nparts):
		if p not in halos: halos[p] = Set()


#def construct_parts_halos(graph, nparts, parts, halos):
#	for node, neighbors in graph.iteritems():
#		mypart = parts[node]
#		if mypart not in halos: halos[mypart] = Set()
#		for neighbor in neighbors:
#			if parts[neighbor] != mypart: halos[mypart].add(neighbor)
#	for p in range(0, nparts):
#		if p not in halos: halos[p] = Set()


def compute_halo_based_weights(nparts, parts, halos, weights):
	## compute num halo cells
	nhalocells = {}
	for part, hcells in halos.iteritems():
		nhalocells[part] = len(hcells)
	for p in range(0, nparts):
		if p not in nhalocells: nhalocells[p] = 0
	print "*** nhalocells: ", nhalocells
	## compute num cells
	ncells = {}
	for node, part in parts.iteritems():
		if part not in ncells: ncells[part] = 0
		ncells[part] += 1
	for p in range(0, nparts):
		if p not in ncells: ncells[p] = 0
	## compute num total cells
	ntotcells = map(add, ncells.values(), nhalocells.values())
	print "*** ncells: ", ncells
	print "*** total cells: ", ntotcells
	print "*** total cells stats: ",
	print_stats(ntotcells)
	pweights = {}
	for part in range(0, nparts):
		#pweights[part] = compute_part_weight(nhalocells[part], ncells[part], nparts)
		pweights[part] = compute_part_weight_model(nhalocells[part], ncells[part], nparts)
	for node, part in parts.iteritems():
		weights[node] = pweights[part]


## equal weights to cells and halo cells
def compute_part_weight(nhalocells, ncells, nparts):
	if ncells != 0: return 1. + float(nhalocells) / ncells
	else: return 1.


## modeled weights to cells and halo cells
def compute_part_weight_model(nhalocells, ncells, nparts):
#	localcompcost = ncells * A / (300 * pow(nparts, 0.84))
#	halocompcost = nhalocells * A / (300 * pow(nparts, 0.84))
#	commcost = nhalocells * B / (300 * pow(nparts, 0.84))
	localcompcost = ncells * A / (pow(nparts, 0.84))
	halocompcost = nhalocells * A / (pow(nparts, 0.84))
	commcost = nhalocells * B / (pow(nparts, 0.84))
	addweight = (commcost + halocompcost) / ncells
	totalweight = (localcompcost / ncells) + addweight
	return totalweight * 10000


def print_stats(arr):
	pdarr = pd.Series(arr)
	asum = pdarr.sum()
	amin = pdarr.min()
	amax = pdarr.max()
	amean = pdarr.mean()
	astd = pdarr.std()
	print 'total: ' + str(asum) + ', min: ' + str(amin) + ', max: ' + str(amax) + ', mean: ' + str(amean) + ', std: ' + str(astd) + ', imbalance: ' + str(1 - (float(amin) / amax))


def write_weighted_hypergraph(filename, data, weights, nnodes, nnets, npins):
	if nnodes != nnets: print 'warning: nnodes and nnets do not match: %d, %d\n' % (nnodes, nnets)
	ff = open(filename, 'w')
	header = '0 %d %d %d 1\n' % (nnodes, nnets, npins)
	ff.write(header)
	for net, pins in sorted(data.iteritems()):
		record = ''
		for p in pins: record += str(p) + ' '
		record += '\n'
		ff.write(record)
	record = ''
	for w in weights: record += str(w) + ' '
	record += '\n'
	ff.write(record)
	ff.close()


def write_parts(filename, parts):
	ff = open(filename, 'w')
	for node, part in sorted(parts.iteritems()):
		record = str(part) + '\n'
		ff.write(record)
	ff.close()


infile, partfile, outfile, nparts, niter = parse_arguments(sys.argv[1:])
weight_factor = lambda x: int(x * 1)

## the first run
data = {}	## mapping { net -> list of nodes in net }
parts = {}	## mapping { node -> partition number }
#graph = {}	## mapping { node -> list of node's neighbors }
## the obtained data and graph always have base index as 0
nnodes, nnets, npins = read_partitioned_hypergraph(infile, partfile, data, parts)
halos = {}
construct_parts_halos(data, nparts, parts, halos)

## compute weights	
weights = map(int, np.ones(nnodes))
compute_halo_based_weights(nparts, parts, halos, weights)
weights = map(weight_factor, weights)
totpweights = {}
calculate_parts_weights(nparts, parts, weights, totpweights)
print "*** total part weights: ", totpweights
tmpoutfile = '.tmp.hypergraph'
write_weighted_hypergraph(tmpoutfile, data, weights, nnodes, nnets, npins)

infile = tmpoutfile
partfile = tmpoutfile + '.part.' + str(nparts)

err = 1e20	## something big
random.seed()
logfile = open('.tmp.log', 'w')

for i in range(0, niter):
	## perform partitioning
	print "============== ", i, " ============="
	seed = random.randint(1, 400000)
	cmd = '~/sw/patoh/build/Linux-x86_64/patoh ' + tmpoutfile + ' ' + str(nparts)
	opts = 'UM=O PQ=S BO=C PA=3 RA=6 SD=' + str(seed)
	cmd += ' ' + opts
	print cmd
	try:
		ret = subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=logfile)
	except subprocess.CalledProcessError:
		print 'error: patoh failed!'
		sys.exit(2)
	else:
		if ret != 0:
			print 'error: something bad happened while running patoh!'
			sys.exit(2)
	parts = {}
	read_partitioning(partfile, parts)

	## compute new halos based on new partitioning
	halos = {}
	construct_parts_halos(data, nparts, parts, halos)

	## compute the new weights
	weights = map(int, np.ones(nnodes))
	compute_halo_based_weights(nparts, parts, halos, weights)
	weights = map(weight_factor, weights)

	totpweights = {}
	calculate_parts_weights(nparts, parts, weights, totpweights)

	newerr = 1. - (float(min(totpweights.values())) / max(totpweights.values()))
	differr = err - newerr
	if differr > 0: accept = True
	else:
		prob = np.exp(differr / TEMPERATURE)
		print '*** Acceptance probability: ' + str(prob)
		if random.random() < prob: accept = True
		else: accept = False;
	
	print "*** ERROR: " + str(newerr) + " [ " + str(err) + ' ]'
	if not accept:
		print '*** rejecting partitioning'
		continue

	err = newerr
	print '*** accepting partitioning'	

	print "*** total part weights: ", totpweights
	print "*** part weight stats: ",
	print_stats(totpweights.values())

	if i == niter: break

	## write new weighted hypergraph
	write_weighted_hypergraph(tmpoutfile, data, weights, nnodes, nnets, npins)

## write out the final partitioning. writing the hypergraph file is not needed ...
totpweights = {}
calculate_parts_weights(nparts, parts, weights, totpweights)
partoutfile = outfile + '.part.' + str(nparts)
write_parts(partoutfile, parts)

logfile.close()
