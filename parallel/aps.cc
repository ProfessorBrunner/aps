/**
 * @file   aps.cc
 * @brief  Calculates Angular Power Spectrum from data files.
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

//#define APS_DEBUG

#include <stdlib.h>

#include <iostream>

#include "OverdensityMap.h"


int main(int argc, char **argv) {
  if (argc != 1){
    std::cout << "Expected 1 argument" << std::endl
              << "Usage:" << std::endl
              << "aps [overdensity.fits] [bandpowers.bin]" << std::endl;
  }
  OverdensityMap map;
  map.LoadFromFile(argv[1]);

  return EXIT_SUCCESS;
}