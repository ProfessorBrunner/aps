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

include ../PATHS

CC=gcc

CFLAGS+=-O3 -fopenmp -m64 -I$(MKLROOT)/include -I$(HEALPIXROOT)/include
LDFLAGS+=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group
LDFLAGS+=-L$(HEALPIXROOT)/lib
LDLIBS+=-ldl -lpthread -lm -lchealpix -lcfitsio

SOURCES=KL_spectrum.c tools.c angular_power_spectrum.c math.c
LIBS=angular_power_spectrum.h

PREPROCESSOR_DEFS=
ifdef APS_KL_COMPRESSION
PREPROCESSOR_DEFS+=-DAPS_KL_COMPRESSION
endif

all: KL_spectrum

KL_spectrum:  $(SOURCES) $(LIBS)
	$(CC) $(PREPROCESSOR_DEFS) $(CFLAGS) $(SOURCES) $(LIBS) -o $@ $(LDFLAGS) $(LDLIBS)

KL_spectrum_output_test: $(SOURCES) $(LIBS)
	$(CC) $(PREPROCESSOR_DEFS) -DAPS_OUTPUT_TEST  $(CFLAGS) $(SOURCES) $(LIBS) -o $@ $(LDFLAGS) $(LDLIBS)
	#-DAPS_SUPPRESS_MATRIX_STDOUT

#Haven't needed to use this
initialize: .FORCE
	$(MKLROOT)/bin/mklvars.sh intel64


#################################################################################
# Profile
#################################################################################
PROFILE_FLAGS=-pg -g -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls
PERF_FLAGS=-lprofiler -ltcmalloc -fno-inline-functions-called-once

profile: .FORCE
	$(CC) $(CFLAGS) $(SOURCES) $(LIBS) -o KL_spectrum $(LDFLAGS) $(LDLIBS) $(PROFILE_FLAGS)

gperf: .FORCE $(FITS) $(DAT)
	$(CC) $(CFLAGS) $(SOURCES) $(LIBS) -o KL_spectrum $(LDFLAGS) $(LDLIBS) $(PERF_FLAGS)

#################################################################################
# Misc
#################################################################################
.FORCE:

clean: .FORCE
	-rm KL_spectrum
	-rm KL_spectrum_output_test
