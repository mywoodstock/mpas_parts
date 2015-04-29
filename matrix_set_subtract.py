##
 #  Project:
 #
 #  File: matrix_set_subtract.py
 #  Created: Mar 09, 2014
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

import sys, getopt, sets

def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, "a:b:o:h")
		if len(opts) < 3:
			raise getopt.GetoptError("Give arguments")
	except getopt.GetoptError:
		print 'matrix_set_substract.py -a <matrix1> -b <matrix2> -o <outputfile>'
		print 'performs set operation (a - b)'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'matrix_set_substract.py -a <matrix1> -b <matrix2> -o <outputfile>'
			sys.exit()
		elif opt in ("-a"):
			inputfile1 = arg
		elif opt in ("-b"):
			inputfile2 = arg
		elif opt in ("-o"):
			outputfile = arg
	return inputfile1, inputfile2, outputfile


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


def matrix_subtract(data1, data2, result, size):
	mat1 = sets.Set(data1)
	mat2 = sets.Set(data2)
	mat3 = mat1.difference(mat2)
	result += list(mat3)

infile1, infile2, outfile = parse_arguments(sys.argv[1:])
data1 = []
data2 = []
size1, nnz1, header1 = read_input_row_col_val(infile1, data1)
size2, nnz2, header2 = read_input_row_col_val(infile2, data2)
print size1, nnz1, size2, nnz2
if size1 != size2:
	print 'error: sizes of the two matrices do not match'
	sys.exit(2)

result = []
matrix_subtract(data1, data2, result, size1)
result = sorted(result, key = lambda entry: entry[1])
result = sorted(result, key = lambda entry: entry[0])

write_result(outfile, size1, result, header1)
