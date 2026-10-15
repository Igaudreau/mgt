#################################
print """

   LEM Info!
   This script is designed to give useful stats for the LEM method given various inputs

"""
##################################

from math import *
# Put this statement at the top of each Python script to import the
# most commonly-used parts of cantera:

ke = 5	# subgrid kinetic energy
npts = 32	# finest LES resolution (grids per half reaction length)
vis = 0.0036	# viscosity
#Nl = 5		# kolmogorov constant (1D)
#Cl = 1.0	# LEM constant (1D)
Ceps = 1.0	# reportedly ranges from 0.2 to 1.0 (controls kolmogorov scale)
Sct = 0.7	# controls stirring time
Cv = 0.01	# also controls stirring time
filterratio = 1	# if the turbulent viscosity filter size is different

dx = 1.0/npts
u_prime = sqrt(2.0*ke/3.0)
Re = u_prime*dx/vis
print "Re = %.4g" % Re

Re = pow(Re,0.75)

eps = (Ceps*pow(ke,1.5))/(dx)
ks = pow((pow(vis,3)/eps),0.25)




#ks = Nl*dx/Re
Nl = ks*Re/dx

#DTurb = u_prime*dx/Cl

#old
#vt = DTurb*Sct
#Cv = vt/(4*dx*sqrt(ke))

#new
vt = Cv*(filterratio*dx*sqrt(ke))
DTurb = vt/Sct
Cl = DTurb/(u_prime*dx)
print "D_Turb = %.3g" % DTurb

a = pow(dx/ks,5/3) - 1
b = 1 - pow(ks/dx,4/3)
lam = 54.0*DTurb*(a/b)/(5.0*pow(dx,3))
dt_stir = 1.0/(dx*lam)
print "dt = %.3g" % dt_stir


a = pow(ks,-2.0/3.0) - pow(dx,-2.0/3.0)
b = pow(ks,-5.0/3.0) - pow(dx,-5.0/3.0)
lavg = 2.5*a/b
nparts = (dx/ks)
print "dx = %.3g" % dx
print "ks (LEM) = %.3g" % ks

print "l_avg = %.3g" % lavg
print "min nparts (dx/ks) = %.3g" % nparts
lavg = lavg/2
nparts = dx/lavg
print "min nparts (lavg/2) = %.3g" % nparts
print ""
print "Cl = %.3g" % Cl
print "Nl = %.3g" % Nl






