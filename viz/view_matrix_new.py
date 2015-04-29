import sys, getopt
import numpy as np
import pandas as pd
import matplotlib as mp
from sets import Set
from matplotlib import pyplot as plt
from mpltools import color

mp.rcParams['axes.linewidth'] = 0.05


mycolor1 = color.color_mapper((0, 5), cmap='Oranges')
mycolor2 = color.color_mapper((0, 5), cmap='Blues')
mycolor3 = color.color_mapper((0, 5), cmap='Greens')
mycolor4 = color.color_mapper((0, 5), cmap='pink')

## colors
myblue		= '#3355A9'
mylblue		= '#7799ff'
myred		= '#ac2318'
mylred		= '#ff675c'
myyellow	= '#d2a41a'
mygreen		= '#8fdf9d'
mylgreen	= '#8fdf9d'
mygray		= '#999999'

usage = 'view_matrix_new.py -i <inputfile> -t <mm|sparse|adj|hyper|patoh> [-s <size>] [-p <partitionfile>] [-v]'

def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, "hi:t:s:p:v")
		if len(opts) < 2:
			raise getopt.GetoptError("Give arguments")
	except getopt.GetoptError:
		print usage
		sys.exit(2)
	view = False
	size_present = False
	partitions = False
	for opt, arg in opts:
		if opt == '-h':
			print usage
			sys.exit()
		elif opt in ("-i"):
			inputfile = arg
		elif opt in ("-t"):
			mtype = arg
		elif opt in ("-s"):
			size = int(arg)
			size_present = True
		elif opt in ('-p'):
			partfile = arg
			partitions = True
		elif opt == '-v':
			view = True
		else:
			print "error: invalid option specified"
			raise getopt.GetoptError("Wrong arguments")
	if not size_present: size = 0
	if not partitions: partfile = None
	return inputfile, mtype, size, partfile, view

def read_matrix(infile, mtype, partfile=None):
	ff = open(infile)
	line = ff.readline()
	skip = 0
	while line.startswith("%"):
		line = ff.readline()
		skip += 1
	words = line.strip().split()
	skip += 1
	if mtype == "sparse" or mtype == "mm":			## sparse/matrix market format
		if len(words) != 3:
			print "error: input file not in sparse or mm format"
			sys.exit(2)
		nx = int(words[0])
		ny = int(words[1])
		nnz = int(words[2])
		sindex = 1
		ff.close()
		mat = pd.read_table(infile, sep='\s+', header=None, comment='%', skiprows=skip, names=['x', 'y', 'v'])
	elif mtype == "hyper" or mtype == "patoh" or mtype == "adj":	## adjacency matrix
		if mtype == "adj":
			if len(words) != 2:
				print "error: input file not in adj format"
				sys.exit(2)
			sindex = 1
			nx = int(words[0])
			ny = nx
			nnz = int(words[1])
		elif mtype == "hyper" or mtype == "patoh":
			words = line.strip().split()
			if len(words) != 4:
				print "error: input file not in hyper or patoh format"
				sys.exit(2)
			sindex = int(words[0])
			nvertices = int(words[1])
			nnets = int(words[2])	## hyperedges
			nx = nvertices
			ny = nnets
			nnz = int(words[3])		## normal edges
		offset = 1 - sindex
		count = 0
		tmat = []
		while True:
			line = ff.readline()
			if not line: break
			count += 1
			words = map(int, line.strip().split())
			for word in words:
				if [count, word+offset, 1] not in tmat:
					tmat.append([count, word+offset, 1])
		ff.close()
		mat = pd.DataFrame(tmat, columns=['x', 'y', 'v'])
	if partfile:
		tmpparts = []
		ff = open(partfile)
		while True:
			line = ff.readline()
			if not line: break
			words = map(int, line.strip().split())
			tmpparts += words
		parts = []
		for i, r in mat.iterrows():
			parts.append(tmpparts[r['x']-1])
#		tparts = pd.read_table(partfile, header=None, names=['part'])
#		parts = []
#		for i, r in mat.iterrows():
#			parts.append(tparts['part'].ix[r['x']-1])
		mat = mat.join(pd.DataFrame(parts, columns=['part']))
	print nx, ny, nnz
	return nx, ny, nnz, mat

def scale_matrix(nx, ny, mat, size):
	if nx == ny:
		sizex = size
		sizey = size
	elif nx < ny:
		sizex = max(int(size * float(size) / ny), 1)
		sizey = size
	else:
		sizex = size
		sizey = max(int(size * float(size) / nx), 1)
	dx = int(np.ceil(float(nx) / sizex))
	dy = int(np.ceil(float(ny) / sizey))
	newmat = {}
	maxx = 0
	maxy = 0
	for i, r in mat.iterrows():
		## bin along x and y
		bx = r['x'] / dx
		by = r['y'] / dy
		maxx = max(maxx, bx)
		maxy = max(maxy, by)
		if (bx, by) not in newmat: newmat[(bx, by)] = r['v']
		else: newmat[(bx, by)] += r['v']
	print maxx, maxy
	mat = []
	for (kx, ky), v in newmat.iteritems():
		mat.append([kx, ky, v])
	return maxx, maxy, pd.DataFrame(mat, columns=['x', 'y', 'v'])

def scale_part_matrix(nx, ny, mat, size):
	if nx != ny:
		print 'error: partitions not valid for non-square matrices'
		sys.exit(2)
	parts = Set()
	for i, r in mat.iterrows(): parts.add(r['part'])
	nparts = len(parts)
	print nx, nparts
	## first divide into nparts bins
	dpx = int(np.ceil(float(nx) / nparts))
	dpy = dpx
	## then subdivide further to get about size bins
	size = int(np.ceil(float(size) / nparts))
	dx = int(np.ceil(float(dpx) / size))
	dy = dx
	pmap = {}
	for i, r in mat.iterrows(): pmap[r['x']] = r['part']
	## bin along x and y
	PBINS = {}
	for i in range(1, nx+1):
		PBINS[i] = int(pmap[i] * size + int(np.ceil(float(i % dpx) / dx) + 1))
	maxx = max(PBINS.values())
	maxy = maxx
	newmat = {}
	tpsz = {}
	partmax = {}
	for i, r in mat.iterrows():
		bx = PBINS[r['x']]
		by = PBINS[r['y']]
		if (bx, by) not in newmat: newmat[(bx, by)] = [r['v'], r['part']]
		else: newmat[(bx, by)][0] += r['v']
		if r['part'] not in tpsz: tpsz[r['part']] = Set()
		tpsz[r['part']].add(bx)
		if r['part'] not in partmax: partmax[r['part']] = 0
		partmax[r['part']] = max(partmax[r['part']], bx)
	mat = []
	for (kx, ky), [v, part] in newmat.iteritems(): mat.append([kx, ky, v, part])
	partsizes = {}
	partsum = 0
	for k, v in tpsz.iteritems():
		partsizes[k] = len(v)
		partsum += len(v)
	return  partsum, partsum, pd.DataFrame(mat, columns=['x', 'y', 'v', 'part']).sort('part'), partsizes

def cluster_matrix(mat, partsz):
	part_offs = {}
	for p, size in sorted(partsz.iteritems()):
		if p == 0: part_offs[p] = 0
		else: part_offs[p] = part_offs[p-1] + partsz[p-1]
	## old node to new node index mapping
	newnodeindex = {}
	partindices = np.zeros(len(part_offs))
	mat = mat.sort('part').reindex()
	for i, r in mat.iterrows():
		if r['x'] not in newnodeindex:
			newnodeindex[r['x']] = partindices[r['part']] + part_offs[r['part']]
			partindices[r['part']] += 1
	for i, r in mat.iterrows():
		r['x'] = newnodeindex[r['x']]+1
		r['y'] = newnodeindex[r['y']]+1
	return mat


def construct_plot(nx, ny, mat, outfile):
	plt.gca().set_aspect('equal', 'box')
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	plt.gca().yaxis.set_major_locator(plt.NullLocator())
	maxv = pd.Series.max(mat['v'])
	ss = mat['v'].apply(float)
	if 'part' in mat: cc = mat['part']
	else: cc = myblue
	plt.scatter(mat['x'], mat['y'], c=cc, cmap='Dark2', alpha=0.7, marker='s', s=ss, lw=0, edgecolor='w')
	plt.xlim([0, nx])
	plt.ylim([ny, 0])
	plt.savefig(outfile, bbox_inches='tight')

def construct_hinton(nx, ny, mat, outfile):
	plt.gca().set_aspect('equal', 'box')
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	plt.gca().yaxis.set_major_locator(plt.NullLocator())
	maxv = max(mat['v'])
	for i, r in mat.iterrows():
		color = 'black'
		size = r['v']
		rect = plt.Rectangle([r['y'], r['x']], size, size, facecolor=color, linewidth=0, edgecolor=color, alpha=float(r['v'])/maxv)
		plt.gca().add_patch(rect)
	plt.xlim([0, nx])
	plt.ylim([ny, 0])
	plt.savefig(outfile, bbox_inches='tight')


infile, mtype, size, partfile, view = parse_arguments(sys.argv[1:])
if size == 0: size = 1000		## size determines max dimension size when not square
nx, ny, nnz, mat = read_matrix(infile, mtype, partfile)
if not partfile:
	sx, sy, mat = scale_matrix(nx, ny, mat, size)
else:
	sx, sy, mat, partsz = scale_part_matrix(nx, ny, mat, size)
	mat = cluster_matrix(mat, partsz)
if partfile: outfile = partfile + '.pdf'
else: outfile = infile + '.pdf'
construct_plot(sx, sy, mat, outfile)
