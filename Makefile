all: Driver.o EvolveProtostellar.o StarProperties.o \
	protostellar_interface.o \
	get_beta_table.o get_derived_values.o \
	interp2dlin.o update_D_mass.o \
	update_luminosity.o update_radius.o \
	update_stage.o
	gfortran Driver.o EvolveProtostellar.o StarProperties.o \
		     protostellar_interface.o \
			 get_beta_table.o get_derived_values.o \
			 interp2dlin.o update_D_mass.o \
			 update_luminosity.o update_radius.o \
			 update_stage.o -o protostellar_evolution

EvolveProtostellar.o : protostellar_interface.o StarProperties.o
	gfortran -c EvolveProtostellar.F90

Driver.o : Driver.F90 EvolveProtostellar.o
	gfortran -c Driver.F90

StarProperties.o : StarProperties.F90
	gfortran -c StarProperties.F90

protostellar_interface.o : protostellar_interface.F90
	gfortran -c protostellar_interface.F90

interp2dlin.o : interp2dlin.F90
	gfortran -c interp2dlin.F90

get_beta_table.o : get_beta_table.F90
	gfortran -c get_beta_table.F90

get_derived_values.o : protostellar_interface.o 
	gfortran -c get_derived_values.F90

update_radius.o : protostellar_interface.o
	gfortran -c update_radius.F90 

update_stage.o : update_stage.F90
	gfortran -c update_stage.F90

update_D_mass.o : update_D_mass.F90
	gfortran -c update_D_mass.F90

update_luminosity.o : update_luminosity.F90
	gfortran -c update_luminosity.F90

clean :
	rm -f *.o *.mod
