# Quadratic Estimation of Angular Power Spectrum
#
# Brett Hayes
# Professor Robert J. Brunner
# Laboratory for Cosmological Data Mining
# University of Illinois Urbana-Champagne
#
# http://lcdm.astro.illinois.edu/code/aps.html
# https://github.com/ProfessorBrunner/aps

# This Makefile is assuming we are compiling on a 64bit system
# This code uses the mkl library for fast operations
# mkl flag maker
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

CC=gcc
MKLROOT=/opt/intel/composer_xe_2013_sp1.2.144/mkl

CFLAGS=-O3 -fopenmp -m64 -I$(MKLROOT)/include
LDFLAGS=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group
#With /opt/intel/composer_xe_2011_sp1.11.339/mkl these flags worked:
#LDFLAGS=-Wl,-Bstatic -Wl,--start-group -L$(MKLROOT)/lib/intel64  -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -liomp5 -Wl,--end-group -Wl,-Bdynamic
LDLIBS=-ldl -lpthread -lm -lchealpix -lcfitsio


SOURCES=KL_spectrum.c tools.c angular_power_spectrum.c math.c
LIBS=angular_power_spectrum.h

all: KL_spectrum

KL_spectrum:  $(SOURCES) $(LIBS)
	$(CC) $(CFLAGS) $(SOURCES) $(LIBS) -o $@ $(LDFLAGS) $(LDLIBS)

initialize: .FORCE
	$(MKLROOT)/bin/mklvars.sh intel64 

#################################################################################
# Misc
#################################################################################
.FORCE:

clean: .FORCE
	rm KL_spectrum