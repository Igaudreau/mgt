Installation.
Install gsl in /usr/local Make sure to compile with -O3
Add  the lbrary to LD_LIBRARY_PATH in .cshrc
setenv LD_LIBRARY_PATH /usr/local/lib:${LD_LIBRARY_PATH}

usr.c commands
pari <value>	set PAR[i]  to <value>
ppari 		print PAR[i]
react netfile	In mg_torch sets up MESA reaction network in netfile.
react iso	In mg_g iso is either iso7 or iso13.
      		In mg_tim/usr13 iso is iso13 or tim13
		In mg_tim/usr21 is tim21
		In mg_tim/usr46 is tim46
		In mg_torch/ options are tim07, tim13, tim21, file.net
		Use noreact first if a network exists
		where file.net is a mesa file
noreact	   	turns off reaction source terms
eos   	eos   	super or ideal (not in mg_torch)
testeos		Compares eos with Nadyo etc		
cons		Calculates the largest abundance from mass fraction
		conservation (default). Switched off by nocons.
cns value       Sets the source term Courant number (default 0.1).
cns1 <value>    sub steps = 100 + int(dt*srate/cns1)
cns2 <value>    controls integration time  in odesolver via
                dt = cn2/srate where srate is the source rate.
                Note must have cns2 < cns1. Only applied for
		ilev < solev.
solev <value>   sets solev.	
dtp value	Sets the print interval for usrmon
tp_old value	Sets the time of last print for usrmon
god    		Use linearised Godunov
hll		Use HLL
kt		Use Kurganov-Tadmor
keps 		switch on turbulence
kmin <value>	minimum turbulent energy
epsmin <value>	minimum turulent dissipation
mtl <value>	maximum turbulent length scale
none            Does not use GSL integrator, just standard MG method
rk2 	   	2nd order Runge Kutta
rk2imp	   	2nd order implicit Runge Kutta
msbdf	   	Poncy variable order variable step implicit etc
nsub  		Sets  number  of  substeps  for GLS  solver  (not
		relevant  for mg_torch)
screen	   	screen nuclear reactions (default on)
noscreen   	turn off screen
revscreen       turns on screening for certain reverse reactions in
                iso7 network (default off for iso7 on for iso13)
doneut 		in mg_torch switch on neutrino losses. nodoneut switches
		off
doweak 		in mg_torch switch on weak rections. nodoweakk switches
		off		
hydrostatic	constant density temperature
dynamic		turns hydrostatic off (default)
isothermal	no temperature udpdate. noisothermal switches off.
print 		prints presures etc with given density and temperature
                at given position. and abundances at that position.
pfcell pos	Prints information in finest grid cell at pos.
pphi		returns phi in steady density equation
shock 		Computes post shock state given right state and
                shock velocity. Abundances taken from given position.
setsteady	sets solution to postshock state for a shock with given
                state
dosteady	calculates ZND solution. Must have a model to get upstream
		state
atmos  		Calculates atmosphere and outputs according to info
		that it asks for
core  		Calculates a WD core and outputs according to info
		that it asks for		
setint name
mode size	sets interpolation from xq files. mode=1 dns, mode = 2
     		steady, mode = 3 atmosphere
		for reading models must do
		read model
		setint
		read model
atminit	filename
		creates a temperature interpolant for an atmosphere using
		filename_ad/iso_rhp/T.xq		
testatminit x	prints atmosphere density at x
shockpos id nx	finds shock position, pressure, and speed by searching
	    	from right boundary in dimension id. nx is other
		coordinate value.
sonpos		Prints position of sonic point for id = 1, nx = 0
smooth nav cmd  Smooths the file cmd.xq over nav  points and writes it
       	   	to cmd_smooth.xq
phaseplan	creates phase<>.xq for given input. Can be used to plot
		phasese portraits		
soncross <name> uses files name_mach_x.xq and files name_u_x.xq to
	 	calculate upstream sound crossing time.
maxxa <name>	Prints max abundance and position in name_xa_<iso>.xq
maldatconv name reads steady.dat file and outputs corresponding .xq
	   	files. Network must have been set up.
malconv name    Converts Malcolm's d-kappa lc-kappa in name.dat to xq
		files dkout.xq and dlcout.xq
timconv name	Converts Timmes name.dat file to nuclear energy release
		vs time in se.xq.
rieman		Riemann problem resolved state		
rat value	Set rat to value
curv value	Set curvature to value
depth value	Set depth to value 
tequi value	Set tequi to value. Used to control temperature in
      		hydrostatic calculations
rho_0 value	Set rho_0 (initial density) to value
T_0 value   	Set T_0 (initial temperature) to value
rho_s value	Set rho_s to value
T_s value   	Set T_s  to value
setmassfrac     Sets mass fraction for the given isotope in
                xl < x < xr. This does not check that mass fractions
		add up to 1.
setrprog pos	Set reaction progress variable to zero for x > pos
setprimv        sets solution to have input rho, u1 or
		temperature in xl < x < xr.
cleansl dir pos Sets the solution to the initial value for x < pos.
       		Used to despike shocks.
cleansr dir pos	Sets the solution to the initial value for x > pos.
vshift value	Shifts the velocity by value.

in main.c
newdom		Now can be used to increase size of domain. New cells
 		use usrqis to get solution at their position.
		Can be used to follow a detonation for a long time. In
		1D only works for one processor
sett <value>    Sets time to value.

plot commands
logab		Uses log abscissa in readxq.
xshift <value>  Shifts abscissa value by <value>.
xsign  		Changes abscissa sign.
plot variables
xa<i> 		plots the mass fraction of isotope i (1..niso).
xan<index>	plots the mass fraction of isotope with name (e.g. he4)
inten		Use internal energy per unit mass for pressure decode..
rri		Reaction coefficient RR[i] for isotope i (1..niso).
sumxa		Plots sum of mass fractions (should be 1.0).
xsoi		mass source for isotope i (1..niso) per unit mass.
		i = niso + 1 gives the energy source per unit mass.
rti		Reaction timescale for isotope i (1..niso).
		i = niso + 1 gives  energy timescale
binden		Binding energy per unit mass.		
abar		Mean atomic mass number.
zbar		Mean atomic number.
usrsum.c	Uses reaction rate for reaction progress variable to
                compute detonation speed.
dtrat 		plots time step reduction factor in usrgls (stored in
                cell).
nsub 		sets min no of sub timesteps (dynamically determined and
		stored in cell can be plotted with by choosing NSUB as
		plot variable.
phi		numerator in steady rho equation. Used for finding
		singular point.

Particles
In evol_part.c  particle velocity  is set to  fluid velocity.  Time at
which particle crosses  shock is put in props[0]  in usrprate.c. Fluid
density,  velocity  and  temperature  put  in  props[1,  props[2]  and
props[3] in usrprate.c

mkpart pos	Puts particle at pos
printpart	prints particle properties
monp		puts particle temperature and density in trho.traj
nomonp		closes trho.traj






%%%%%%%%%%%%%%%%%%%
ZND example

$read timmesmix_34
mon znd
curv 1e-08 dosteady 0.9e9 1.1e9 1.0e9 verbose/grot bisect
			Dlow Dhigh Dguess

atmos
base density, min density, temperature, stellar mass, radius, shell mass, adiabatic/isothermal, 
output filename 

atmos
mix_atmos_x_13.dat
mix_atmos_xy_13.dat




