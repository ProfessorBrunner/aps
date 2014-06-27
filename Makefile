# Quadratic Estimation of Angular Power Spectrum
#
# Brett Hayes
# Professor Robert J. Brunner
# Laboratory for Cosmological Data Mining
# University of Illinois Urbana-Champagne
#
# http://lcdm.astro.illinois.edu/code/aps.html
# https://github.com/ProfessorBrunner/aps


all: shared_memory distributed_memory

shared_memory: .FORCE
	$(MAKE) -C ./shared_memory

distributed_memory: .FORCE
	$(MAKE) -C distributed_memory

#################################################################################
# Test
#################################################################################

TEST_NSIDE=8
FITS=data/$(TEST_NSIDE)_53918_lcdm.fits
DAT=data/CL_$(TEST_NSIDE)_lcdm.dat

test: profile  $(FITS) $(DAT)
	./shared_memory/KL_spectrum $(FITS) $(DAT)

%.fits: data test/create_mock.py

%.dat: data test/create_mock.py

data: .FORCE
	./test/create_mock.py ./test/catalog.dat ./data $(TEST_NSIDE)

profile: .FORCE
	$(MAKE) -C ./shared_memory profile


#################################################################################
# Misc
#################################################################################
.FORCE:

clean: .FORCE
	$(MAKE) -C ./shared_memory clean
	$(MAKE) -C ./distributed_memory clean

clean_data: .FORCE
	rm data/*