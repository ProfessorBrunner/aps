# Quadratic Estimation of Angular Power Spectrum
#
# Brett Hayes
# Professor Robert J. Brunner
# Laboratory for Cosmological Data Mining
# University of Illinois Urbana-Champagne
#
# http://lcdm.astro.illinois.edu/code/aps.html
# https://github.com/ProfessorBrunner/aps


include PATHS

all: shared_memory distributed_memory

shared_memory: .FORCE
	$(MAKE) -C ./shared_memory

distributed_memory: .FORCE
	$(MAKE) -C distributed_memory

#################################################################################
# Test
#################################################################################

TEST_NSIDE=8
NUM_PROC=1
#FITS=$(TEST_NSIDE)_53918_lcdm
#BANDS=CL_$(TEST_NSIDE)_lcdm
FITS=32_node_compareb_1000000000
BANDS=32_node_compareb_1000000000
DISTRIBUTED_TEST_DIR=data/test_distributed_$(BANDS)
SHARED_TEST_DIR=data/test_shared_$(BANDS)

FITS_PATH=data/$(FITS).fits
BANDS_PATH=data/$(BANDS).bands

COMPARE_FLAGS=-o data/kl_graphs-2/ --heat-plot

test_shared: shared_memory/KL_spectrum_output_test $(FITS_PATH) $(BANDS_PATH)
	rm -rf $(SHARED_TEST_DIR)
	./shared_memory/KL_spectrum_output_test $(FITS_PATH) $(BANDS_PATH)
	./test/cat_signal.bash $(SHARED_TEST_DIR)
	$(PYTHON) ./test/compare_test_directories.py data/standard $(SHARED_TEST_DIR) $(COMPARE_FLAGS)

test_distributed: distributed_memory/aps_test $(FITS_PATH) $(BANDS_PATH)
	rm -rf $(DISTRIBUTED_TEST_DIR)
	mpirun -n $(NUM_PROC) ./distributed_memory/aps_test $(FITS_PATH) $(BANDS_PATH)
	./test/cat_signal.bash $(DISTRIBUTED_TEST_DIR)
	$(PYTHON) ./test/compare_test_directories.py data/standard $(DISTRIBUTED_TEST_DIR) $(COMPARE_FLAGS)

test_debug: shared_memory/KL_spectrum_output_test  distributed_memory/aps_test $(FITS_PATH) $(BANDS_PATH)
	@echo
	@echo "Running tests and creating output"
	@echo "Clearing test output directories"
	-rm -rf $(DISTRIBUTED_TEST_DIR)
	-rm -rf $(SHARED_TEST_DIR)

	./shared_memory/KL_spectrum_output_test $(FITS_PATH) $(BANDS_PATH)
	mpirun -n $(NUM_PROC) ./distributed_memory/aps_test $(FITS_PATH) $(BANDS_PATH)

	# ./test/cat_signal.bash $(SHARED_TEST_DIR)
	# ./test/cat_signal.bash $(DISTRIBUTED_TEST_DIR)
	$(PYTHON) ./test/compare_test_directories.py $(SHARED_TEST_DIR) $(DISTRIBUTED_TEST_DIR)  $(COMPARE_FLAGS)

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

standard: .FORCE
	@echo
	@echo "Copying into standard"
	@echo
	rm -rf ./data/standard
	cp -R $(SHARED_TEST_DIR) data/standard


%.fits: data test/generate_inputs.py
	@:

%.bands: data test/generate_inputs.py
	@:

data: .FORCE
	$(PYTHON) ./test/generate_inputs.py -c ./test/catalog.dat ./data $(TEST_NSIDE)

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
	-rm -rf data
