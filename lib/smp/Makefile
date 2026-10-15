# smp/Makefile
#
.SILENT:
# SMP components
SMP_OBJ= \
smp.o 

SMP_LIB=smp.a
# Link smp
smp_always:	 $(SMP_OBJ)
	echo "Creating smp.a library $(Comp_opts) $(C_opts)"
	rm -f $(SMP_LIB)
	ar cr $(SMP_LIB) $(SMP_OBJ)
	rm *.o

#COMPILE
smp.o: smp.c; $(C_comp) $(Comp_opts) $(C_opts) -c smp.c
