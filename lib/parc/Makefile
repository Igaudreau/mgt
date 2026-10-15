# parc/Makefile
.SILENT:
# Libraries

# locations
parc_obj=.
compile=.

#OBJECTS
OBJ= \
$(parc_obj)/parc.o \
$(parc_obj)/dam.o  

# Link parc
all:	$(MAIN) $(OBJ)
	echo "Creating parc$(MP).a library $(Comp_opts) $(C_opts)"
	rm -f parc$(MP).a
	ar cr parc$(MP).a parc.o dam.o
	rm *.o

# MPI
mpi:
	make all MP="mpi" MP_CC="$(mpi_CC) $(Comp_opts) $(C_opts) \
		-DPARC_MPI" \
		MP_LD="$(mpi_CC)" 

# SMP
smp:
	make all MP="smp" MP_CC="$(C_comp) $(Comp_opts) $(C_opts) \
		-DPARC_SMP" \
		MP_LD="$(C_comp)" 

both:
	make mpi
	make smp

#COMPILE
$(parc_obj)/parc.o: parc.c;  $(MP_CC) -c parc.c
$(parc_obj)/dam.o: dam.c; $(MP_CC) -c dam.c
