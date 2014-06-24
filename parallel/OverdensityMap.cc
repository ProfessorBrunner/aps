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

#include <iostream>

#include "chealpix.h"

#include "OverdensityMap.h"

OverdensityMap::OverdensityMap() {
  nside_ = 0;
  bins_ = 0;
  total_galaxies_ = 0;
  omega_ = 0;
}

OverdensityMap::~OverdensityMap() {
  free(healpix_map_);
  free(ra_);
  free(dec_);
}

void OverdensityMap::LoadFromFile(char *file) {
  std::cout << "Loading overdensity from: " << file << std::endl;
  char ordering[10], coords[1];
  long nside;
  healpix_map_ = read_healpix_map(file, &nside, coords, ordering);


  assert(coords[0] == 'C');
  assert(strcmp(ordering, "NESTED") == 0);
  assert(nside > 0);
  nside_ = nside;
}