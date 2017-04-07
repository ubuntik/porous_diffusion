#!/opt/local/bin/python2.7
#

def is_number(s):
	try:
		float(s) # for int, long and float
	except ValueError:
		return False
	return True

n = 200
inf = open('data.vtk', 'r+')
outf = open('plot.dat', 'w+')

data = list()
for i in range(0, n):
	data.append(str(str(i) + " "))

cnt = 0
lst = list(inf)
for l in lst:
	if is_number(l):
		data[cnt] = '%s %f' % (data[cnt], float(l))
		cnt = cnt + 1
		if cnt == n:
			cnt = 0

for l in data:
	outf.write('%s\n' % l)

