import sys, getopt

def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, "hi:t:s:v")
		if len(opts) < 2:
			raise getopt.GetoptError("Give arguments")
	except getopt.GetoptError:
		print 'view_matrix.py -i <inputfile> -t <type:imat/gmat/hmat> -s <size>'
		sys.exit(2)
	for opt, arg in opts:
		view = False
		size_present = False
		if opt == '-h':
			print 'view_matrix.py -i <inputfile> -t <type:imat/gmat/hmat> -s <size>'
			sys.exit()
		elif opt in ("-i"):
			inputfile = arg
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
	return inputfile, mtype, size, view

def my_open_matrix(infile, mat, mtype, symmetrize):
	ff = open(infile)
	line = ff.readline()
	while line.startswith("%"): line = ff.readline()
	if mtype == "gmat":			## graph.info style
		words = line.strip().split()
		nnodes = int(words[0])
		nedges = int(words[1])
		count = 1
		while True:
			line = ff.readline()
			if not line: break
			words = line.strip().split()
			if count not in mat: mat[count] = {}
			for word in words:
				mat[count][int(word)] = 1.0
				if symmetrize:
					if int(word) not in mat: mat[int(word)] = {}
					mat[int(word)][count] = 1.0
			count += 1
	elif mtype == "imat":
		words = line.strip().split()
		if len(words) == 3:
			nnodes = int(words[0])
			nedges = int(words[2])
			while True:
				line = ff.readline()
				if not line: break
				words = line.strip().split()
				row = int(words[0])
				col = int(words[1])
				val = float(words[2])
				if row not in mat: mat[row] = {}
				mat[row][col] = val
				if symmetrize:
					if col not in mat: mat[col] = {}
					mat[col][row] = val
	ff.close()
	return nnodes, nedges


infile, mtype, size, view = parse_arguments(sys.argv[1:])

from summon import matrix

if mtype == "hmat":
	mat = matrix.Matrix()
	matrix.open_matrix(infile, mat, format = mtype, symmetrize = True)
else:
	mat = {}
	nnodes, nedges = my_open_matrix(infile, mat, mtype, True)

mini = minj = maxi = maxj = 1
for i, cols in mat.iteritems():
	maxi = max(maxi, i)
	mini = min(mini, i)	
	for j, val in cols.iteritems():
		maxj = max(maxj, j)
		minj = min(minj, j)
maxy = maxx = max(maxi, maxj)

print mini, minj, maxi, maxj
if size == 0: size = maxx

black = (0, 0, 0)
white = (255, 255, 255)

outfile = open(infile + ".svg", "w")
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
for i, cols in mat.iteritems():
	for j, val in cols.iteritems():
		#x = i
		#y = j
		x = (float(i - mini + 1) / (maxi - mini + 1)) * size - 1
		y = (float(j - minj + 1) / (maxj - minj + 1)) * size - 1
		#elem = "<circle cx='%d' cy='%d' r='0.5' stroke-opacity='1.000000' " % (x, y)
		elem = "<rect x='%d' y='%d' height='1' width='1' " % (x, y)
		elem += "stroke='rgb(%d,%d,%d)' stroke-width='0' />\n" % black		## points
		outfile.write(elem)
footer = "</g></g></g>\n</g></g></g>\n</svg>"
outfile.write(footer)
outfile.close()

if view:
	viewer = matrix.MatrixViewer(mat, winsize = (size, size), title = "Sparse Matrix", bgcolor=(1, 1, 1))
	viewer.show()
