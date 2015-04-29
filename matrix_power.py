##
 #  Project:
 #
 #  File: matrix_power.py
 #  Created: Mar 09, 2014
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

import sys, getopt

def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, "i:o:p:")
		if len(opts) < 3:
			raise getopt.GetoptError("Give arguments")
	except getopt.GetoptError:
		print 'matrix_power.py -i <inputfile> -o <outputfile> -p <power>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'matrix_power.py -i <inputfile> -o <outputfile> -p <power>'
			sys.exit()
		elif opt in ("-i"):
			inputfile = arg
		elif opt in ("-o"):
			outputfile = arg
		elif opt in ("-p"):
			power = int(arg)
	return inputfile, outputfile, power


def read_input_row_col_val(filename, data):
	ff = open(filename)
	head = ff.readline()
	header = ""
	while head.startswith("%"):
		header += head
		head = ff.readline()
	headwords = head.strip().split()
	if len(headwords) == 2:
		size = int(headwords[0])
		nnz = int(headwords[1])
	elif len(headwords) == 3:
		size = int(headwords[0])
		size2 = int(headwords[1])
		if size != size2:
			print "error: input matrix is not square"
			sys.exit(2)
		nnz = int(headwords[2])
	count = 0
	while True:
		line = ff.readline()
		if not line: break
		line = line.strip()
		count += 1
		words = line.split()
		data += [(int(words[0]), int(words[1]), int(words[2]))]	# assuming vals == 1 only
	if count != nnz:
		print "error: count != nnz"
		exit
	ff.close()
	return size, nnz, header


def write_result(filename, size, result, header = ""):
	ff = open(filename, 'w')
	head = "%s%d\t%d\t%d\n" % (header, size, size, len(result))
	ff.write(head)
	for i in range(0, len(result)):
		record = "%d\t%d\t%d\n" % (result[i][0], result[i][1], result[i][2])
		ff.write(record)
	ff.close()


def matrix_square(data, result, size):
	# since matrix is symmetric, col[i] == row[i]
	# there is not symmetric data included. avoiding its construction.
	row_data = sorted(data, key = lambda entry: entry[0])
	col_data = sorted(data, key = lambda entry: entry[1])
	rrcount = 0
	rccount = 0
	for i in range(1, size + 1):
		rrow = []
		rcol = []
		while rrcount < len(row_data) and row_data[rrcount][0] == i:
			rcol += [row_data[rrcount][1]]
			rrcount += 1
		while rccount < len(col_data) and col_data[rccount][1] == i:
			rrow += [col_data[rccount][0]]
			rccount += 1
		row = sorted(rcol + rrow)
		crcount = 0
		cccount = 0
		for j in range(1, i):
			crow = []
			ccol = []
			while crcount < len(row_data) and row_data[crcount][0] == j:
				ccol += [row_data[crcount][1]]
				crcount += 1
			while cccount < len(col_data) and col_data[cccount][1] == j:
				crow += [col_data[cccount][0]]
				cccount += 1
			col = sorted(ccol + crow)
			rowindex = 0
			colindex = 0
			index = 0
			val = 0
			while rowindex < len(row) and colindex < len(col):
				if(row[rowindex] == col[colindex]):
					val += 1
					rowindex += 1
					colindex += 1
				elif(row[rowindex] < col[colindex]): rowindex += 1
				elif(row[rowindex] > col[colindex]): colindex += 1
			if(val > 0): result += [(j, i, 1)]


def matrix_cube(data, result, size):
	row_data = sorted(data, key = lambda entry: entry[0])
	col_data = sorted(data, key = lambda entry: entry[1])
	sq = []
	matrix_square(data, sq, size)
	sqrow_data = sorted(sq, key = lambda entry: entry[0])
	sqcol_data = sorted(sq, key = lambda entry: entry[1])
	sqrcount = 0
	sqccount = 0
	for i in range(1, size + 1):
		sqrow_r = []
		sqrow_c = []
		while sqrcount < len(sqrow_data) and sqrow_data[sqrcount][0] == i:
			sqrow_c += [sqrow_data[sqrcount][1]]
			sqrcount += 1
		while sqccount < len(sqcol_data) and sqcol_data[sqccount][1] == i:
			sqrow_r += [sqcol_data[sqccount][0]]
			sqccount += 1
		row = sorted(sqrow_c + sqrow_r)
		if len(row) == 0: continue
		crcount = 0
		cccount = 0
		for j in range(1, i):
			col_r = []
			col_c = []
			while crcount < len(row_data) and row_data[crcount][0] == j:
				col_c += [row_data[crcount][1]]
				crcount += 1
			while cccount < len(col_data) and col_data[cccount][1] == j:
				col_r += [col_data[cccount][0]]
				cccount += 1
			col = sorted(col_c + col_r)
			if len(col) == 0: continue
			rowindex = 0
			colindex = 0
			val = 0
			while rowindex < len(row) and colindex < len(col):
				if(row[rowindex] == col[colindex]):
					val += 1
					rowindex += 1
					colindex += 1
				elif(row[rowindex] < col[colindex]): rowindex += 1
				elif(row[rowindex] > col[colindex]): colindex += 1
			if(val > 0): result += [(j, i, 1)]


inpfile, outfile, power = parse_arguments(sys.argv[1:])
data = []
size, nnz, header = read_input_row_col_val(inpfile, data)
print size, nnz
result = []
if power == 2:
	matrix_square(data, result, size)
elif power == 3:
	matrix_cube(data, result, size)
else:
	print "only powers 2 and 3 are recognized for now"
	sys.exit(2)
result = sorted(result, key = lambda entry: entry[0])
write_result(outfile, size, result, header)
