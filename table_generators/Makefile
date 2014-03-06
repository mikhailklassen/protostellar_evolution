# DESCRIPTION:
#
# This makefile is used to build all the tools of the polytropes_table directory
#
# These tools are used to:
#       + generate MODELDATA:      a table listing polytropic parameters for
#                                  different polytropic stars
#
#       + generate BetaTable.dat:  a table listing in three columns the mass,
#                                  volume-averaged beta, and derivate of log(beta)
#                                  with respect to log(m). It lists these in
#                                  segments, with each segment corresponding a
#                                  polytropic index n. They are computed for the
#                                  same values of n listed in MODELDATA and in the
#                                  same order.
#
# CALLING SEQUENCE:
#
#       make modeldata             Creates the executable to generate MODELDATA
#
#       make betatable             Creates the executable to generate BetaTable.dat
#
#       make clean                 Deletes the .o and .mod files
#
# WRITTEN: Mikhail Klassen         (McMaster University)
#          15 / July / 2010
#
# EMAIL: klassm@mcmaster.ca
#

# Module macros:
#
MODULES = astro_constants.f90 \
          polytools.f90
        # Modules that contain procedures and must therefore
        # be compiled together with project
CMODS   = polytools.f90
M_OBJS  = make_modeldata.f90 \
          solveLE.f90
B_OBJS  = make_beta_table.f90 \
          solve_beta.f90

# Definitions
#
RM := rm -rf

# Configuration:
#
FC = gfortran
CFLAGS = -O3

# The projects to be built
#
all : modeldata betatable

modeldata : modules modeldata_OBJS

betatable : modules betatable_OBJS

betatable_OBJS : $(CMODS) $(B_OBJS)
	@echo 'Compiling beta table generator...'
	$(FC) $(CFLAGS) $(CMODS) $(B_OBJS) -o make_betatable
	@echo 'Success.'

modeldata_OBJS : $(M_OBJS)
	@echo 'Compiling modeldata table generator...'
	$(FC) $(CFLAGS) $(M_OBJS) -o make_modeldata
	@echo 'Success.'

modules : $(MODULES)
	@echo 'Compiling modules:'
	@echo '   $(MODULES)'
	$(FC) $(CFLAGS) -c $(MODULES)
	@echo 'Done.'

interptest : $(CMODS) examples/get_beta_table.f90 interp2dlin.f90
	@echo 'Compiling beta tables interpolator test'
	@echo '... get_beta_table.f90 from the examples folder'
	@echo '... inter2dlin.f90' 
	$(FC) $(CFLAGS) $(CMODS) examples/get_beta_table.f90 interp2dlin.f90 -o interptest
	@echo 'Done.'

# Delete .mod and .o files
clean : 
	-$(RM) *.mod *.o
    
# Also delete executables
realclean :
	-$(RM) *.mod *.o make_betatable make_modeldata
