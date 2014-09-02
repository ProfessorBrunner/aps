# Quadratic Estimation of Angular Power Spectrum
# Parallel implementation

include /usr/local/conf/ElemVars
CC=gcc
CXX=g++
MKLROOT=/opt/intel/composer_xe_2013_sp1.2.144/mkl
HEALPIXROOT=/usr/local/src/Healpix_3.11/

CFLAGS+=-I$(HEALPIXROOT)/include
LDFLAGS+=-L$(HEALPIXROOT)/lib
LDLIBS+=-ldl -lm -lcfitsio -lchealpix

SOURCES=AngularPowerSpectrum.cc aps.cc OverdensityMap.cc BandPower.cc
LIBS=AngularPowerSpectrum.h BandPower.h OverdensityMap.h

PREPROCESSOR_DEFS=
ifdef APS_KL_COMPRESSION
PREPROCESSOR_DEFS+=-DAPS_KL_COMPRESSION
endif

all: aps

aps:  $(SOURCES) $(LIBS)
	$(CXX) $(PREPROCESSOR_DEFS) $(CFLAGS) ${ELEM_COMPILE_FLAGS} $(SOURCES) $(LIBS) -o $@ $(LDFLAGS) ${ELEM_LINK_FLAGS} $(LDLIBS) ${ELEM_LIBS}

aps_test: $(SOURCES) $(LIBS)
	$(CXX) $(PREPROCESSOR_DEFS) -DAPS_OUTPUT_TEST $(CFLAGS) ${ELEM_COMPILE_FLAGS} $(SOURCES) $(LIBS) -o $@ $(LDFLAGS) ${ELEM_LINK_FLAGS} $(LDLIBS) ${ELEM_LIBS}

#Haven't needed to use this
initialize: .FORCE
	$(MKLROOT)/bin/mklvars.sh intel64

#################################################################################
# Misc
#################################################################################
.FORCE:

clean: .FORCE
	-rm aps aps_test