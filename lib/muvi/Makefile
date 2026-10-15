# Make MUVI (LINUX version)
#
.SILENT:
OBJ = $(GEN) $(DEV)

GEN =\
chr.o \
fchr.o \
cont.o \
dr3.o \
dr3mat.o \
dr3solid.o \
dr3wire.o \
drp.o \
ds.o \
csg.o \
model.o \
modelv.o \
cuser.o \
nam.o \
sys.o \
vis.o 

DEV=\
ps.o \
x11.o 

# Link muvi
muvi:	$(OBJ)
	echo "Creating muvi.a library $(Comp_opts) $(C_opts)" 
	touch muvi.a
	rm muvi.a
	ar cr muvi.a $(OBJ)
	rm *.o

# GEN Compilation
chr.o: chr.c; $(C_comp) -c $(Comp_opts) $(C_opts) chr.c
fchr.o: fchr.c; $(C_comp) -c $(Comp_opts) $(C_opts) fchr.c 
cont.o: cont.c; $(C_comp) -c $(Comp_opts) $(C_opts) cont.c 
dr3.o: dr3.c; $(C_comp) -c $(Comp_opts) $(C_opts) dr3.c 
dr3mat.o: dr3mat.c; $(C_comp) -c $(Comp_opts) $(C_opts) dr3mat.c 
dr3solid.o: dr3solid.c; $(C_comp) -c $(Comp_opts) $(C_opts) dr3solid.c 
dr3wire.o: dr3wire.c; $(C_comp) -c $(Comp_opts) $(C_opts) dr3wire.c 
drp.o: drp.c; $(C_comp) -c $(Comp_opts) $(C_opts) drp.c 
ds.o: ds.c; $(C_comp) -c $(Comp_opts) $(C_opts) ds.c 
csg.o: csg.c; $(C_comp) -c $(Comp_opts) $(C_opts) csg.c 
model.o: model.c; $(C_comp) -c $(Comp_opts) $(C_opts) model.c 
modelv.o: modelv.c; $(C_comp) -c $(Comp_opts) $(C_opts) modelv.c 
cuser.o: cuser.c; $(C_comp) -c $(Comp_opts) $(C_opts) cuser.c 
nam.o: nam.c; $(C_comp) -c $(Comp_opts) $(C_opts) nam.c 
sys.o: sys.c; $(C_comp) -c $(Comp_opts) $(C_opts) sys.c 
vis.o: vis.c; $(C_comp) -c $(Comp_opts) $(C_opts) vis.c 

# DEV Comp
ps.o: ps.c; $(C_comp) -c $(Comp_opts) $(C_opts) ps.c 
x11.o: x11.c; $(C_comp) -c $(Comp_opts) $(C_opts) x11.c 
