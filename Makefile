#*********************************************************************************************************
#**                                     Developped by Tomas Tamfal 					**
#**                                       Version 0.1 April 2017	                               	**
#*********************************************************************************************************


#**********************************************************************************************************
#            Makefile_start: Compiles the code with different compilers and different parameters
#**********************************************************************************************************
SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .F90 .o

#********** Compiler and options choice **********
#FC = pgf90
FC = gfortran

#------- GFORTRAN ------
FFLAGS = -O3 
#FFLAGS = -O0 -Og -Wall -Wextra -pedantic -fimplicit-none -fcheck=all  -fbacktrace  
#FFLAGS = -O0 -Wunused-variable
#*********** Library linking ******************************
LDLIBS = #

#********** Choice of simulation parameters **********
#DENSITY = ALPHABETAGAMMA
DENSITY = HERNQUIST
#DENSITY = SIS
#DENSITY = NFW
#DENSITY = TEST

INTEGRATION = TRAPEZ

TIMESTEPS = FIXED
#TIMESTEPS = VARIABLE

MESH = LIN_MESH
#MESH = LOG_MESH


#********** Object files **********
modules = parameters.o variables.o

objects = main.o mesh.o integration.o density.o angularmomentum.o \
	dynamicalfrictionforce.o lnlambda.o velocitydistribution.o  

target = ../stopping_time.run


#********** Commands **********
all: def $(target) cl

$(target): $(modules) $(objects)
	$(FC) $(FFLAGS) -o $@ $^ ${LDLIBS} 
	
%.o: %.F90
	$(FC) $(FFLAGS) -c  $<  ${LDLIBS}
	
%.o: choice.def

.PHONY: cl clean

cl:
	rm -f *~ *.o *.mod OUT choice.def
	
clean:
	rm -f *~ *.o *.mod *data fort* OUT
		
def:
	echo \#define $(DENSITY) >> choice.def
	echo \#define $(INTEGRATION) >> choice.def	
	echo \#define $(TIMESTEPS) >> choice.def
	echo \#define $(MESH) >> choice.def









