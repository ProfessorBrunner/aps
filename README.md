# Quadratic Estimation of Angular Power Spectrum

* * *

Professor Robert J. Brunner
Laboratory for Cosmological Data Mining
University of Illinois Urbana-Champagne

### Original shared memory implementation
Brett Hayes
http://lcdm.astro.illinois.edu/code/aps.html
https://github.com/ProfessorBrunner/aps

### Distributed memory implementation 
Work done during:
#### The Passionate On Parallel REU
at The University of Illinois Urbana-Champaign
by Joy Hill and Alex Warren
Under the advisement of Matias Carrasco Kind

* * *

# Installation
Code was developed with Ubuntu 12.04 and greater
but should be compatible with other Unix enviroments

## Dependencies

#### Shared memory implementation

[CFITSIO 3.370](http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
[HEALPix 3.11](http://healpix.jpl.nasa.gov/)
[MKL](https://software.intel.com/en-us/intel-mkl) or another library for the
BLAS/LAPACK routines.

#### Shared memory implementation

[CFITSIO 3.370](http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
[HEALPix 3.11](http://healpix.jpl.nasa.gov/)
[Elemental 0.84](http://libelemental.org/releases/0.84/index.html)

## Configuration
#### Shared memory implementation

After installing the dependencies:
	cd shared_memory
	cp Makefile.template Makefile
and edit the Makefile to find all the libraries.