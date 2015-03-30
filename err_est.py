#!/usr/bin/python

import sys
from optparse import OptionParser
from itertools import izip
from os import listdir
from os.path import isfile, join
from math import sqrt, fabs

max_val = 0

def isfloat(d):
	try:
		float(d)
	except ValueError:
		return False
	return True

def filter_lines(lines, mod):
	olines = []
	i = 0
	for line in lines:
		i = i + 1
		if not isfloat(line):
			olines.append(line)
			i = i - 1
		elif (i - 1) % mod == 0:
			olines.append(float(line))
	return olines

def C_norm(n1, n2):
	res = fabs(n1*n1 - n2*n2)
	return format(res / fabs(max(n1, n2, key=abs)), '.6f')

def L2_norm(n1, n2):
	res = fabs(n1*n1 - n2*n2)
	dev = sqrt(n1*n1 + n2*n2) / 2.0
	return format(res / dev, '.6f')

def exec_counting(indirs, outdir):
	files = []
	for cdir in indirs:
		files.append([f for f in listdir(cdir) if isfile(join(cdir,f))])
#	num_files = 1
	num_files = len(files[0])
	for i in range(0,num_files):
		file1 = open(indirs[0]+files[0][i], 'r')
		file2 = open(indirs[1]+files[1][i], 'r')
		ofile = open(outdir+files[0][i], 'w+')

		lines1 = filter_lines(file1.readlines(), 10)
		lines2 = filter_lines(file2.readlines(), 1)

		print len(lines1), len(lines2)

		for line1, line2 in izip(lines1, lines2):
			if not isfloat(line2):
				ofile.write(line2)
			elif line1 == line2 == 0:
				ofile.write("0.000000\n")
			else:
				ofile.write(str(C_norm(line1, line2))+"\n")

def parse_args():
	parser = OptionParser()
	parser.add_option("-i", "--input", dest="indir",
                  nargs = 2, help="get data from dirs")
	parser.add_option("-o", "--output", dest="outdir",
                  help="write result to dir")

	(options, args) = parser.parse_args()
	return options.indir, options.outdir

def main():
	indirs, outdir = parse_args()
	exec_counting(indirs, outdir)

if __name__ == "__main__":
	main()
