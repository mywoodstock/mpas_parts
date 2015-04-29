import sys, getopt, random

def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, "hi:p:t:s:v")
		if len(opts) < 3:
			raise getopt.GetoptError("Give arguments")
	except getopt.GetoptError:
		print 'view_matrix_parts.py -i <inputfile> -p <partfile> -t <type:imat/gmat/hmat> -s <size>'
		sys.exit(2)
	for opt, arg in opts:
		view = False
		size_present = False
		if opt == '-h':
			print 'view_matrix_parts.py -i <inputfile> -p <partfile> -t <type:imat/gmat/hmat> -s <size>'
			sys.exit()
		elif opt in ("-i"):
			inputfile = arg
		elif opt in ("-p"):
			partfile = arg
		elif opt in ("-t"):
			mtype = arg
		elif opt in ("-s"):
			size = int(arg)
			size_present = True
		elif opt == '-v':
			view = True
		else:
			print "error: invalid option specified"
			raise getopt.GetoptError("Wrong arguments")
	if not size_present: size = 0
	return inputfile, partfile, mtype, size, view

def read_matrix(infile, partfile, mat, mtype):
	pp = open(partfile)
	parts = []
	for partline in pp: parts += map(int, partline.strip().split())
	pp.close()
	ff = open(infile)
	line = ff.readline()
	while line.startswith("%"): line = ff.readline()
	words = line.strip().split()
	if mtype == "gmat":
		nnodes = int(words[0])
		nedges = int(words[1])
		count = 1
		nparts = 0
		while True:
			line = ff.readline()
			if not line: break
			#partline = pp.readline()
			#if not partline: break
			#part = int(partline.strip().split()[0])
			part = parts[count - 1]
			if part not in mat:
				mat[part] = {}
				nparts += 1
			words = line.strip().split()
			if count not in mat[part]: mat[part][count] = {}
			for word in words:
				mat[part][count][int(word)] = 1.0
			count += 1
	elif mtype == "hmat":
		nnodes = int(words[1])
		nedges = int(words[3])
		count = 1
		nparts = 0
		while True:
			line = ff.readline()
			if not line: break
			part = parts[count - 1]
			if part not in mat:
				mat[part] = {}
				nparts += 1
			words = line.strip().split()
			if count not in mat[part]: mat[part][count] = {}
			for word in words:
				mat[part][count][int(word)] = 1.0
			count += 1
	else:
		print "error: unsupported file format"
		ff.close()
		#pp.close()
		sys.exit(2)
	ff.close()
	return nnodes, nedges, nparts

def cluster_matrix(mat, nparts, partsizes, newmat):
	sorted_partsizes = sorted(partsizes.items())
	part_offsets = {}
	for p, size in sorted_partsizes:
		if p == 0: part_offsets[p] = 0
		else: part_offsets[p] = part_offsets[p - 1] + sorted_partsizes[p - 1][1]
	# create old to new node index mapping
	newnodeindex = {}
	tempmat = {}
	for p, part in mat.iteritems():
		tempmat[p] = sorted(part.items())
		i = 0
		for n, neighbors in tempmat[p]:
			i += 1
			newnodeindex[n] = i + part_offsets[p]
	for p, plist in tempmat.iteritems():
		newmat[p] = {}
		for n, neighbors in plist:
			newn = newnodeindex[n]
			newmat[p][newn] = {}
			for k, v in neighbors.iteritems():
				newmat[p][newn][newnodeindex[k]] = v


infile, partfile, mtype, size, view = parse_arguments(sys.argv[1:])
mat = {}
if mtype == "gmat" or mtype == "hmat":
	nnodes, nedges, nparts = read_matrix(infile, partfile, mat, mtype)
else:
	print "only 'gmat' and 'hmat' mtype supported for now"
	sys.exit()

mini = minj = maxi = maxj = 1
partsizes = {}
for p, part in mat.iteritems():
	if p not in partsizes: partsizes[p] = len(part)
	for i, cols in part.iteritems():
		maxi = max(maxi, i)
		mini = min(mini, i)	
		for j, val in cols.iteritems():
			maxj = max(maxj, j)
			minj = min(minj, j)
maxy = maxx = max(maxi, maxj)

newmat = {}
cluster_matrix(mat, nparts, partsizes, newmat)

if size == 0: size = maxx
# construct colors
black = (0, 0, 0)
white = (255, 255, 255)
partcolors = {}
for i in range(0, nparts):
	partcolors[i] = (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))

outfile = open(partfile + ".svg", "w")
header = "<?xml version='1.0' encoding='UTF-8'?>\n"
header += "<!DOCTYPE svg PUBLIC '-//W3C//DTD SVG 1.1//EN' "
header += "'http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd'>\n"
header += "<svg width='%d' height='%d' xmlns='http://www.w3.org/2000/svg' version='1.1'>\n" % (size, size)
header += "<g style='font-family: courier'>\n"
header += "<g transform='translate(0, 0)'>\n"
header += "<g transform='scale(1, 1)'>\n"
header += "<rect x='0' y='0' width='%d' height='%d' " % (size, size)
header += "fill='rgb(%d,%d,%d)'/>\n" % white		## background
header += "<g transform='scale(1.000000, 1.000000)'>\n"
header += "<g transform='translate(0, 0)'>\n"
header += "<g stroke-width='1.000000'>\n"
outfile.write(header)
for p, part in newmat.iteritems():
	for i, cols in part.iteritems():
		for j, val in cols.iteritems():
			x = (float(i - mini + 1) / (maxi - mini + 1)) * size - 1
			y = (float(j - minj + 1) / (maxj - minj + 1)) * size - 1
			elem = "<rect x='%d' y='%d' height='1' width='1' " % (x, y)
			elem += "fill='rgb(%d,%d,%d)' stroke-width='0' />\n" % partcolors[p]		## points
			outfile.write(elem)
footer = "</g></g></g>\n</g></g></g>\n</svg>"
outfile.write(footer)
outfile.close()

if view:
	viewer = matrix.MatrixViewer(mat, winsize = (size, size), title = "Sparse Matrix", bgcolor=(1, 1, 1))
	viewer.show()
