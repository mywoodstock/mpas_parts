##
 #  Project:
 #
 #  File: generate_unweighted_partitions.py
 #  Created: June 25, 2014
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

import sys, getopt
import subprocess
import random


def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, 'i:n:')
		if len(opts) != 2:
			raise getopt.GetoptError('Give arguments')
	except getopt.GetoptError:
		print 'compute_unweighted_partitions.py -i <hypergraphfile> -n <nparts>'
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-i'): inputfile = arg
		elif opt in ('-n'): nparts = int(arg)
	return inputfile, nparts


infile, nparts = parse_arguments(sys.argv[1:])

random.seed()
seed = random.randint(1, 400000)
cmd = '~/sw/patoh/build/Linux-x86_64/patoh ' + infile + ' ' + str(nparts)
opts = 'UM=O PQ=S BO=C PA=3 RA=6 SD=' + str(seed)
cmd += ' ' + opts
print cmd
try:
	ret = subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError:
	print 'error: patoh failed!'
	sys.exit(2)
else:
	if ret != 0:
		print 'error: something bad happened while running patoh!'
		sys.exit(2)
