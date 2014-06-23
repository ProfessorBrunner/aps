# Quadratic Estimation of Angular Power Spectrum
#
# Brett Hayes
# Professor Robert J. Brunner
# Laboratory for Cosmological Data Mining
# University of Illinois Urbana-Champagne
#
# http://lcdm.astro.illinois.edu/code/aps.html
# https://github.com/ProfessorBrunner/aps


all: sequential

sequential: .FORCE
	$(MAKE) -C ./sequential

parallel: .FORCE
	$(MAKE) -C parallel

#################################################################################
# Test
#################################################################################

TEST_NSIDE=8
FITS=data/$(TEST_NSIDE)_53918_lcdm.fits
DAT=data/CL_$(TEST_NSIDE)_lcdm.dat

test: profile  $(FITS) $(DAT)
	./sequential/KL_spectrum $(FITS) $(DAT)

%.fits: data test/create_mock.py

%.dat: data test/create_mock.py

data: .FORCE
	./test/create_mock.py ./test/catalog.dat ./data $(TEST_NSIDE)

profile: .FORCE
	$(MAKE) -C ./sequential profile


#################################################################################
# Misc
#################################################################################
.FORCE:

clean: .FORCE
	$(MAKE) -C ./sequential clean

clean_data: .FORCE
	rm data/*