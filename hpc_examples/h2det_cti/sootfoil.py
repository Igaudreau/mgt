import csv
from array import *
from math import *
from numpy import *

# First parse command line arguments to get number of steps and marches and export options

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--f", dest="f", default=1)
parser.add_option("--l", dest="l", default=200)
(options, args) = parser.parse_args()

firstFile = int(options.f)
lastFile = int(options.l)


# open file and output the following 1D info
filenum = str(firstFile).zfill(4)
currentFile = './xq/sf' + filenum + '.xq'
data = loadtxt(currentFile,skiprows=3)
head = open(currentFile)

# get number of lines in list
count = sum(1 for line in data)
# need to create arrays this size!
maxsf = [0.0] * count

# loop through all sootfoil files and find maximum for each data point
for j in range(firstFile, lastFile+1):
	filenum = str(j).zfill(4)
	currentFile = './xq/sf' + filenum + '.xq'
	data = loadtxt(currentFile,skiprows=3)

	# get maximum
	for i in range(0, count):
		if (data[i][2] > maxsf[i]):
			maxsf[i] = data[i][2]



# now create new XQ file

outFile = 'sootfoil.xq'
#print num
f = open(outFile,'w')
f.write(head.readline())
f.write(head.readline())
f.write(head.readline())
head.close()
#f.write("export\n")
#f.write("X Y SF\n ")
for i in range(0, count):
	f.write("%s" %data[i][0])
	f.write(" ")
	f.write("%s" %data[i][1])
	f.write(" ")
	f.write("%s" %maxsf[i])
	f.write("\n")
f.close()






