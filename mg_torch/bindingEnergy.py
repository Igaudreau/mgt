from math import *
from numpy import *

# First parse command line arguments to get number of steps and marches and export options

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--e", dest="e", default=0)
parser.add_option("--n", dest="n", default=2)
parser.add_option("--p", dest="p", default=2)
parser.add_option("--m", dest="m", default=4.0026)
(options, args) = parser.parse_args()


e = float(options.e)
n = float(options.n)
p = float(options.p)
mass = float(options.m)

me = 0.000549
mp = 1.007825
mn = 1.008665
A = n+p
BE = 931.5*(e*me + n*mn + p*mp - mass)/A

print "Binding Energy %.4g MeV/nucleon" % BE

BE = 1.6e-6*BE

print "Binding Energy %.4g ergs/nucleon" % BE


BE = 6.02214e23*BE

print "Binding Energy %.4g ergs/g" % BE

