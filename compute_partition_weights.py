##
 #  Project:
 #
 #  File: generate_weighted_patoh.py
 #  Created: Jul 30, 2014
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


TEMPERATURE = 1e-3
A = 2. / 3.

def F(p):
	return np.log(4. * p)



def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, 'i:p:o:y:n:')
		if len(opts) < 5 or len(opts) > 5:
			raise getopt.GetoptError('Give arguments')
	except getopt.GetoptError:
		print 'compute_partitions_weights.py -i <hypergraphfile> -p <partitionfile> -o <outputfile> -y <nhalolayers> -n <nparts>'
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-i'): inputfile = arg
		elif opt in ('-o'): outputfile = arg
		elif opt in ('-p'): partfile = arg
		elif opt in ('-n'): nparts = int(arg)
		elif opt in ('-y'): nhalos = int(arg)
	return inputfile, partfile, outputfile, nparts, nhalos


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
	read_partitioning(partfile, parts)

	return nnodes, nnets, npins


def get_halo_candidates(data, cell, candidates, nlayers):
	if nlayers < 1: return
	for pin in data[cell]:
		candidates.add(pin)
		get_halo_candidates(data, pin, candidates, nlayers - 1)


def construct_parts_halos(data, nparts, parts, nlayers, halos):
	print nparts, ':: constructing halos with', nlayers, 'layers ...'
	## collect cells in partitions
	cells = {}
	for p in range(0, nparts):
		cells[p] = []
	for node, part in parts.iteritems():
		cells[part] += [ node ]
	## get halo candidates for all parts
	for p in range(0, nparts):
		candidates = set()
		for cell in cells[p]:
			get_halo_candidates(data, cell, candidates, nlayers)
		halos[p] = list(candidates - set(cells[p]))


def compute_part_neighbors(halos, nparts, parts, pneighbors):
	for p in range(0, nparts):
		pneighbors[p] = set()
		for hcell in halos[p]:
			pneighbors[p].add(parts[hcell])
			if p == parts[hcell]:
				print 'error: something really bad happened!'
				sys.exit(2)


def compute_model_computation_cost(ncells, nhcells, nparts, local_comps, halo_comps):
	a = A
	f = F(nparts)
	for p in range(0, nparts):
		local_comps[p] = ncells[p] / f
		halo_comps[p] = a * nhcells[p] / f


def compute_model_communication_cost(nhcells, nneighbors, tot_comps, nparts, tot_comms):
	f = F(nparts)
	max_comp = max(tot_comps.values())
	for p in range(0, nparts):
		tot_comms[p] = (nhcells[p] / (nneighbors[p] * f)) + (max_comp - tot_comps[p])


def calculate_parts_total_weights(nparts, parts, weights, pweights):
	for p in range(0, nparts):
		pweights[p] = 0
	for node in range(0, len(weights)):
		part = parts[node]
		pweights[part] += weights[node]


def compute_weights(nnodes, nparts, parts, halos, pneighbors, outfile):
	## compute num halo cells
	nhalocells = {}
	for p in range(0, nparts):
		nhalocells[p] = 0
	for p, hcells in halos.iteritems():
		nhalocells[p] = len(hcells)
	## compute num cells
	ncells = {}
	for p in range(0, nparts):
		ncells[p] = 0
	for node, part in parts.iteritems():
		ncells[part] += 1
	## compute num total cells
	ntotcells = map(add, ncells.values(), nhalocells.values())
	## compute num neighbors
	nneighbors = {}
	for p in range(0, nparts):
		nneighbors[p] = 0
	for p, n in pneighbors.iteritems():
		nneighbors[p] = len(n)
	## compute cost due to computations
	local_comps = {}
	halo_comps = {}
	compute_model_computation_cost(ncells, nhalocells, nparts, local_comps, halo_comps)
	tot_comps = {}
	for p in range(0, nparts):
		tot_comps[p] = local_comps[p] + halo_comps[p]
	## compute cost due to communications
	tot_comms = {}
	compute_model_communication_cost(nhalocells, nneighbors, tot_comps, nparts, tot_comms)
	## calculate the weights to add to each partition cell due to halo cells
	part_add_weights = {}
	part_local_weights = {}
	for p in range(0, nparts):
		part_add_weights[p] = float(halo_comps[p] + tot_comms[p]) / ncells[p]
		part_local_weights[p] = float(local_comps[p]) / ncells[p]	## averaging
	weights = map(int, np.ones(nnodes))
	for node, part in parts.iteritems():
		weights[node] = int((part_local_weights[part] + part_add_weights[part]) * 10000)
		#weights[node] = part_add_weights[part] * 10000
	pweights = {}
	calculate_parts_total_weights(nparts, parts, weights, pweights)
	write_parts_stats(outfile, nparts, ncells, nhalocells, ntotcells, nneighbors, local_comps, halo_comps, tot_comps, tot_comms, pweights)


def write_parts_stats(outfile, nparts, ncells, nhalocells, ntotcells, nneighbors, local_comps, halo_comps, tot_comps, tot_comms, pweights):
	ff = open(outfile, 'w')
	header = 'partition	ncells	nhcells	ntotcells	nneighbors	localcomp	halocomp	totcomp	totcomm	modeledweight\n'
	ff.write(header)
	for p in range(0, nparts):
		record = '%d	%d	%d	%d	%d	%f	%f	%f	%f	%d\n' % (p, ncells[p], nhalocells[p], ntotcells[p], nneighbors[p], local_comps[p], halo_comps[p], tot_comps[p], tot_comms[p], pweights[p])
		ff.write(record)
	ff.close()


## the main part
## for a given partitioning, output:
##	{ partition ncells nhcells ntotcells nneighbors localcomp halocomp totcomp comm modeledweight }

infile, partfile, outfile, nparts, nlayers = parse_arguments(sys.argv[1:])

data = {}	## mapping { net -> list of nodes in net }
parts = {}	## mapping { node -> partition number }
## the obtained data and graph always have base index as 0
nnodes, nnets, npins = read_partitioned_hypergraph(infile, partfile, data, parts)
halos = {}
construct_parts_halos(data, nparts, parts, nlayers, halos)
pneighbors = {}
compute_part_neighbors(halos, nparts, parts, pneighbors)
## compute weights	
compute_weights(nnodes, nparts, parts, halos, pneighbors, outfile)
