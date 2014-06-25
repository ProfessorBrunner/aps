/**
 * @file   OverdensityMap.cc
 * @brief  Used to load overdensity map and header data from a .fits file.
 *
 * @date   June, 2014
 * @Author Robert Brunner
 * @Author Brett Hayes
 * @Author Matias Carrasco-Kind
 * @Author Joy Hill
 * @Author Alex Warren
 *
 * This algorithm was originally developed by Brett Hayes in:
 * The Laboratory for Cosmological Data Mining
 * Professor Robert J. Brunner
 * http://lcdm.astro.illinois.edu/
 * University of Illinois Urbana-Champagne
 *
 * Originally published: http://lcdm.astro.illinois.edu/papers/sdssdr7-aps.html
 *
 * Work continued through the
 * Passionate on Parallel Research Experience for Undergraduates (REU)
 * Orginal code was developed for shared memory machines
 * During the REU, Alex and Joy converted the code to work with distributed
 * memory using MPI.
 */


#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

#include <iostream>

#include "chealpix.h"
#include "fitsio.h"

#include "OverdensityMap.h"

OverdensityMap::OverdensityMap() {
  nside_ = 0;
  bins_ = 0;
  total_galaxies_ = 0;
  omega_ = 0;
  healpix_map_ = NULL;
  ra_ = NULL;
  dec_ = NULL;
}

OverdensityMap::~OverdensityMap() {
  if (healpix_map_) free(healpix_map_);
  if (ra_) free(ra_);
  if (dec_) free(dec_);
}

void OverdensityMap::LoadFromFile(char *file) {
  std::cout << "Loading overdensity from: " << file << std::endl;
  char ordering[10], coords[1];
  long nside;
  healpix_map_ = read_healpix_map(file, &nside, coords, ordering);

  assert(coords[0] == 'C');
  assert(strcmp(ordering, "NESTED") == 0);
  
  LoadFitsKeys(file);

  if (nside_ != nside) {
    std::cout << "read_healpix_map outputs " << nside << "." << std::endl
              << "fits_read_key_lng outputs " <<  nside_ << "." << std::endl
              << "Using " << nside_ << "." << std::endl
              << "This bug may only exist on Alex's system." << std::endl;
  }
  assert(nside_ > 0);
  assert(total_galaxies_ > 0);
}

void OverdensityMap::LoadFitsKeys(char *file) {
  fitsfile *fptr;
  int status=0, hdutype;
  long long_value;

  fits_open_file(&fptr, file, READONLY, &status);
  if (status) fits_report_error(stderr, status);

  fits_movabs_hdu(fptr, 2, &hdutype, &status);
  if (status) fits_report_error(stderr, status);
  assert(hdutype==BINARY_TBL);

  fits_read_key_lng(fptr, "NSIDE", &long_value, NULL, &status);
  nside_ = (int) long_value;
  if (status) {
    std::cerr << "Error reading key: NSIDE" << std::endl;
    fits_report_error(stderr, status);
  } 

  fits_read_key_lng(fptr, "NGALAXY", &long_value, NULL, &status);
  total_galaxies_ = (double) long_value;
  if (status) {
    std::cerr << "Error reading key: NGALAXY" << std::endl;
    fits_report_error(stderr, status);
  }

  fits_close_file(fptr, &status);
  if (status) {
    std::cerr << "Error closing file" << std::endl;
    fits_report_error(stderr, status);
  }

  std::cout << "NSIDE:   " << nside_ << std::endl
            << "NGALAXY: " << total_galaxies_ << std::endl;

  assert(status == 0); //Correctly loaded header data
}