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
 * University of Illinois Urbana-Champaign
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
#include <stdlib.h>

#include <iostream>
#include <cstring>
#include <string>

#include "chealpix.h"
#include "fitsio.h"

#include "OverdensityMap.h"

OverdensityMap::OverdensityMap()
    :  nside_(0),
       bins_(0),
       total_galaxies_(0),
       omega_(0),
       healpix_map_(NULL),
       overdensity_(NULL),
       ra_(NULL),
       dec_(NULL) {}

OverdensityMap::~OverdensityMap() {
  if (healpix_map_) free(healpix_map_);
  if (overdensity_) free(overdensity_);
  if (ra_) free(ra_);
  if (dec_) free(dec_);
}

void OverdensityMap::LoadFromFile(char *file_path) {
  std::cout << std::string(80, '-') << std::endl;
  std::cout << "Loading overdensity from: " << file_path << std::endl;
  char ordering[10], coords[1];
  long nside;

  std::cout << "Getting raw pixel data..." << std::endl;
  healpix_map_ = read_healpix_map(file_path, &nside, coords, ordering);

  assert(healpix_map_ != NULL);
  assert(coords[0] == 'C');
  assert(std::strcmp(ordering, "NESTED") == 0);
  
  LoadFitsKeys(file_path);

  ReadHealpixMap();

  if (nside_ != nside) {
    std::cout << "Mismatched value for nside:" << std::endl
              << "read_healpix_map outputs " << nside << "." << std::endl
              << "fits_read_key_lng outputs " <<  nside_ << "." << std::endl
              << "Using " << nside_ << "." << std::endl
              << "This bug may only exist on Alex's system." << std::endl;
  }
  assert(nside_ > 0);
  assert(total_galaxies_ > 0);
  std::cout << "Finished loading overdensity map." << std::endl;
  std::cout << std::string(80, '-') << std::endl;
}

void OverdensityMap::LoadFitsKeys(char *file_path) {
  std::cout << "Getting FITS keys..." << std::endl;
  fitsfile *fptr;
  int status=0, hdutype;
  long long_value;

  fits_open_file(&fptr, file_path, READONLY, &status);
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

void OverdensityMap::ReadHealpixMap() {
  std::cout << "Reading Healpix map..." << std::endl;
  int total_pixels = nside2npix(nside_);
  double theta, phi;
  bins_ = 0;


  // TODO(Alex): these loops could be parallelized with openmp
  for (int i = 0; i < total_pixels; ++i) {
    if (healpix_map_[i] >= -1.0) ++bins_;
  }

  ra_ = (double *) malloc(bins_*sizeof(double));
  dec_ = (double *) malloc(bins_*sizeof(double));
  overdensity_ = (double *)  malloc(bins_*sizeof(double));

  int j = 0;
  for (int i = 0; i < total_pixels; ++i) {
    if (healpix_map_[i] >= -1.0) {
      pix2ang_nest(nside_, i, &theta, &phi);
      ra_[j] = phi*kRadianToDegree;
      dec_[j] = 90.0 - theta*kRadianToDegree;
      if (ra_[j] < 0.0) ra_[j] += 360.0;
      overdensity_[j] = healpix_map_[i];

      // if (j < 20){
      //   std::cout << ra_[j] << " " << dec_[j] << std::endl;
      // } 
      ++j;
    }
  }
  
  omega_ = bins_ * kSquareDegreePerSphere / total_pixels;

  std::cout << "Loaded " << bins_ << " bins out of "
            << total_pixels << " possible pixels with total area of "
            << omega_ << " square degrees." << std::endl;
}