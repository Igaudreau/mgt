.SILENT:
# GNU
CC=$(CPP_comp)
F90=$(F77_comp)
MP_CC = $(mpi_CC)
OPTS="$(Comp_opts) $(CPP_opts)"
#F90_OPTS = "-O3 --std=legacy "
#F90_OPTS = "-g --std=legacy "
F90_OPTS="$(Fortran_opts)"

# mg program
OBJ=\
linux/main.o \
linux/start.o \
linux/march.o \
linux/mg.o \
linux/mug.o \
linux/adi.o \
linux/grid.o \
linux/utils.o \
linux/para.o \
linux/plot.o \
linux/readxq.o \
linux/mvcmd.o \
linux/modelrw.o \
linux/help.o \
linux/iop.o \
linux/usr.o \
linux/torch.o \
linux/tim13.o \
linux/tim21.o
#linux/fpe.o

#Muvi Library
LIB =$(mg)/lib/muvi.a

# Equilibrium-c package
EQLIB = $(mg)/lib/libeqc.a

# Link mg
all:	$(OBJ)
	touch linux/mg
	rm linux/mg
	echo "Linking mg (new) "
	$(MP_LD) -o linux/mg $(OBJ) \
	$(MP_LIB) $(LIB) $(EQLIB) -lgfortran \
	$(X11_libs) $(system)

mpi:
	make MP="mpi" MP_LIB="$(mg)/lib/parcmpi.a" MP_LD=$(MP_CC)

smp:
	make MP="smp" MP_LIB="$(mg)/lib/parcsmp.a \
	$(mg)/lib/smp.a " \
	MP_LD=$(CC)

#Compile
linux/main.o: main.c; linux/ac $(CC) main $(OPTS)
linux/start.o: start.c; linux/ac $(CC) start $(OPTS)
linux/march.o: march.c; linux/ac $(CC) march $(OPTS)
linux/mg.o: mg.c; linux/ac $(CC) mg $(OPTS)
linux/mug.o: mug.c; linux/ac $(CC) mug $(OPTS)
linux/adi.o: adi.c; linux/ac $(CC) adi $(OPTS)
linux/grid.o: grid.c; linux/ac $(CC) grid $(OPTS)
linux/utils.o: utils.c; linux/ac $(CC) utils $(OPTS)
linux/para.o: para.c; linux/ac $(CC) para $(OPTS)
linux/plot.o: plot.c; linux/ac $(CC) plot $(OPTS)
linux/readxq.o: readxq.c; linux/ac $(CC) readxq $(OPTS)
linux/mvcmd.o: mvcmd.c; linux/ac $(CC) mvcmd $(OPTS)
linux/modelrw.o: modelrw.c; linux/ac $(CC) modelrw $(OPTS)
linux/help.o: help.c; linux/ac $(CC) help $(OPTS)
linux/iop.o: iop.c; linux/ac $(CC) iop $(OPTS)
linux/usr.o: usr.c; linux/ac $(CC) usr $(OPTS)
#linux/fpe.o: fpe.c; linux/ac $(CC) fpe $(OPTS)
linux/torch.o: torch.f90; linux/af $(F90) torch $(F90_OPTS)
linux/tim13.o: tim13.f90; linux/af $(F90) tim13 $(F90_OPTS)
linux/tim21.o: tim21.f90; linux/af $(F90) tim21 $(F90_OPTS)

