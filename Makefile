# Definitions
#
RM := rm -rf

# Configuration:
#
FC = gfortran
CFLAGS = -O3

all: Driver.o EvolveProtostellar.o \
	protostellar_interface.o \
	get_beta_table.o get_derived_values.o \
	interp2dlin.o update_D_mass.o \
	update_luminosity.o update_radius.o \
	update_stage.o
	@echo 'Compiling protostellar evolution code...'
	$(FC) $(CFLAGS) Driver.o EvolveProtostellar.o \
		     protostellar_interface.o \
			 get_beta_table.o get_derived_values.o \
			 interp2dlin.o update_D_mass.o \
			 update_luminosity.o update_radius.o \
			 update_stage.o -o protostellar_evolution
	@echo 'Success.'

EvolveProtostellar.o : protostellar_interface.o 
	$(FC) $(CFLAGS) -c EvolveProtostellar.F90

Driver.o : Driver.F90 EvolveProtostellar.o
	$(FC) $(CFLAGS) -c Driver.F90

protostellar_interface.o : protostellar_interface.F90
	$(FC) $(CFLAGS) -c protostellar_interface.F90

interp2dlin.o : interp2dlin.F90
	$(FC) $(CFLAGS) -c interp2dlin.F90

get_beta_table.o : get_beta_table.F90
	$(FC) $(CFLAGS) -c get_beta_table.F90

get_derived_values.o : protostellar_interface.o 
	$(FC) $(CFLAGS) -c get_derived_values.F90

update_radius.o : protostellar_interface.o
	$(FC) $(CFLAGS) -c update_radius.F90 

update_stage.o : update_stage.F90
	$(FC) $(CFLAGS) -c update_stage.F90

update_D_mass.o : update_D_mass.F90
	$(FC) $(CFLAGS) -c update_D_mass.F90

update_luminosity.o : update_luminosity.F90
	$(FC) $(CFLAGS) -c update_luminosity.F90

clean :
	@echo 'Deleting .o and .mod files'
	$(RM) *.o *.mod

realclean :
	@echo 'Deleting .o and .mod files, and the binary executable'
	$(RM) *.o *.mod protostellar_evolution
