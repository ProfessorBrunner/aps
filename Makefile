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
DAT=data/CL_$(TEST_NSIDE)_lcdm.bands
NUM_PROC=1

test_shared: shared_memory/KL_spectrum_output_test $(FITS) $(DAT)
	./shared_memory/KL_spectrum_output_test $(FITS) $(DAT)
	./test/compare_test_directories.py data/standard data/test_shared_CL_$(TEST_NSIDE)_lcdm

test_distributed: distributed_memory/aps_test $(FITS) $(DAT)
	rm -rf data/test_distributed_CL_$(TEST_NSIDE)_lcdm
	mpirun -n $(NUM_PROC) ./distributed_memory/aps_test $(FITS) $(DAT)
	./test/compare_test_directories.py data/standard data/test_distributed_CL_$(TEST_NSIDE)_lcdm

test_debug: shared_memory/KL_spectrum_output_test  distributed_memory/aps_test $(FITS) $(DAT)
	@echo
	@echo "Running tests and creating output"
	@echo "Clearing test output directories"
	rm -rf data/test_distributed_CL_$(TEST_NSIDE)_lcdm
	rm -rf data/test_shared_CL_$(TEST_NSIDE)_lcdm
	./shared_memory/KL_spectrum_output_test $(FITS) $(DAT)
	mpirun -n $(NUM_PROC) ./distributed_memory/aps_test $(FITS) $(DAT)
	./test/compare_test_directories.py data/test_shared_CL_$(TEST_NSIDE)_lcdm data/test_distributed_CL_$(TEST_NSIDE)_lcdm

distributed_memory/aps_test: .FORCE
	@echo
	@echo "Making distributed angular power spectrum binary"
	@echo
	$(MAKE) -C ./distributed_memory aps_test

shared_memory/KL_spectrum_output_test: .FORCE
	@echo
	@echo "Making shared memory angular power spectrum binary"
	@echo
	$(MAKE) -C ./shared_memory KL_spectrum_output_test

%.fits: data test/generate_inputs.py
	@:

%.bands: data test/generate_inputs.py
	@:

data: .FORCE
	./test/generate_inputs.py -c ./test/catalog.dat ./data $(TEST_NSIDE)

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