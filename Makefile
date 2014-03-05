all: Driver.o EvolveProtostellar.o StarProperties.o \
	MSFits.o polytools.o \
	get_beta_table.o get_derived_values.o \
	interp2dlin.o update_D_mass.o \
	update_luminosity.o update_radius.o \
	update_stage.o
	gfortran Driver.o EvolveProtostellar.o StarProperties.o \
		     MSFits.o polytools.o \
			 get_beta_table.o get_derived_values.o \
			 interp2dlin.o update_D_mass.o \
			 update_luminosity.o update_radius.o \
			 update_stage.o -o protostellar_evolution

EvolveProtostellar.o : MSFits.o StarProperties.o
	gfortran -c EvolveProtostellar.F90

Driver.o : Driver.F90 EvolveProtostellar.o
	gfortran -c Driver.F90

StarProperties.o : StarProperties.F90
	gfortran -c StarProperties.F90

MSFits.o : MSFits.F90
	gfortran -c MSFits.F90

polytools.o : polytools.F90
	gfortran -c polytools.F90

interp2dlin.o : interp2dlin.F90
	gfortran -c interp2dlin.F90

get_beta_table.o : get_beta_table.F90
	gfortran -c get_beta_table.F90

get_derived_values.o : MSFits.o polytools.o
	gfortran -c get_derived_values.F90

update_radius.o : MSFits.o
	gfortran -c update_radius.F90 

update_stage.o : update_stage.F90
	gfortran -c update_stage.F90

update_stage.F90 : polytools.o MSFits.o 
	gfortran -c polytools.F90 

update_D_mass.o : update_D_mass.F90
	gfortran -c update_D_mass.F90

update_luminosity.o : update_luminosity.F90
	gfortran -c update_luminosity.F90

clean :
	rm -f *.o *.mod
