# Quadratic Estimation of Angular Power Spectrum

* * *

Professor Robert J. Brunner  
Laboratory for Cosmological Data Mining  
University of Illinois Urbana-Champaign  

#### Original shared memory implementation
Brett Hayes  
http://lcdm.astro.illinois.edu/code/aps.html  
https://github.com/ProfessorBrunner/aps  

#### Distributed memory implementation
Work done during __The Passionate On Parallel REU__  
at __The University of Illinois Urbana-Champaign__  

By Joy Hill and Alex Warren  
under the advisement of Matias Carrasco Kind  

* * *

# Installation
Code was developed with Ubuntu 12.04 and greater
but should be compatible with other Unix enviroments

## Dependencies

#### Shared memory implementation

 * [CFITSIO 3.370](http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
 * [HEALPix 3.11](http://healpix.jpl.nasa.gov/)
 * [MKL](https://software.intel.com/en-us/intel-mkl) or another library for the
BLAS/LAPACK routines.

#### Distributed memory implementation

 * [CFITSIO 3.370](http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
 * [HEALPix 3.11](http://healpix.jpl.nasa.gov/)
 * [Elemental 0.84](http://libelemental.org/releases/0.84/index.html)

### Input Generation and Testing Scripts
Python scripts are used for creating test data and comparing the results.
The following python packages are used:
 * [numpy](http://www.numpy.org/)
 * [pyfits](http://www.stsci.edu/institute/software_hardware/pyfits)
 * [healpy](http://healpy.readthedocs.org/en/latest/)
 * [matplotlib](http://matplotlib.org/)


 _pyfits_ and _healpy_ are available through pip.

## Configuration
#### Shared memory implementation

After installing the dependencies:

    $ cp PATHS.template PATHS

This file is included in the Makefiles for the two implementations,
edit it to find your libraries.

* * *

# Running
#### Testing implementations against each other
Make sure you you have execute permissions for the the scripts in test/

    $ make test_debug

It may

#### Shared Memory

#### Distributed Memory
