/**
 * @file   AngularPowerSpectrum.cc
 * @brief  Calculates angular power spectrum from an overdensity map of astronomical data.
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
#include <vector>

#include <iostream>
#include <cstring>
#include <string>

#include "elemental-lite.hpp"
#include "chealpix.h"
#include "fitsio.h"

#include "AngularPowerSpectrum.h"


using namespace elem;

//AngularPowerSpectrum::AngularPowerSpectrum() {}

AngularPowerSpectrum::AngularPowerSpectrum(int bins, int bands, 
    double total_galaxies, double omega, Grid &grid)
    :  bins_(bins),
       bands_(bands),
       total_galaxies_(total_galaxies),
       omega_(omega),
       signal_(new std::vector<DistMatrix<double>>(bands)),
       grid_(&grid),
       c_(new double[bands]),
       c_start_(new int[bands]),
       c_end_(new int[bands]),
       ra_(new double[bins]),
       dec_(new double[bins]),
       grid_height_(grid.Height()),
       grid_width_(grid.Width()),
       grid_row_(grid.Row()),
       grid_col_(grid.Col()),
       local_height_(Length(bins*bins, grid_row_, grid_height_)),
       local_width_(Length(bins*bins, grid_col_, grid_width_)) {
}

AngularPowerSpectrum::~AngularPowerSpectrum() {}


void AngularPowerSpectrum::run() {
  CalculateSignal();
}

/**
 * Build Covariance Matrix
 */
void AngularPowerSpectrum::CalculateCovariance() {}


void AngularPowerSpectrum::CalculateSignal() {
  std::vector<double> cos_values(bins_*bins_, 0.0f);
  std::vector<double> previous_previous(bins_*bins_, 1.0f);
  std::vector<double> previous = cos_values;
  std::vector<double> local_sum = cos_values;
  int k = 0;
  for (int ell = 0; ell < c_end_[bands_-1]; ++ell) {
    if (ell == 0) {
      std::vector<double> local_signal(bins_*bins_, 0.0f);

      continue;
    }
    if (ell ==1){
      continue;
    }
    if (ell>=c_end_[k]){
      ell = c_start_[++k];      
      std::vector<double> local_signal(bins_*bins_, 0.0f);
    }
    double coefficient = (2 * ell + 1)/(2 * ell * (ell + 1));


  }
}

/**
 * Karhunen-Loeve Compression
 */
void AngularPowerSpectrum::KLCompression() {}

/**
 * Estimates the APS given Signal & Covariance Matrices
 */
void AngularPowerSpectrum::EstimateC() {}

/**
 * Build Expected Covariance Matrix based on C_l
 * Previously named calculate_difference 
 */
void AngularPowerSpectrum::ExpectedCovariance() {}

/**
 * Build Fisher Matrix & Weighted Average
 */
void AngularPowerSpectrum::CalculateFisher() {}

/**
 * Recalculate C_l from the Fisher Matrix
 */
void AngularPowerSpectrum::RecalculateC_L() {}