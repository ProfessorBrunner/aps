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
#include <math.h>

#include <vector>
#include <iostream>
#include <cstring>
#include <string>
#include <algorithm>
#include <functional>

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
       local_height_(Length(bins, grid_row_, grid_height_)),
       local_width_(Length(bins, grid_col_, grid_width_)) {
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
  bool i_have_it = false;

  signal_ = std::vector<DistMatrix<double>>(bands_, DistMatrix<double>(*grid_));


  std::vector<double> cos_values(local_height_ * local_width_, 0.0f);
  for( Int jLoc = 0; jLoc < local_width_; ++jLoc ) {
      // Form global column index from local column index
      const Int j = grid_col_ + jLoc*grid_width_;
      for( Int iLoc = 0; iLoc < local_height_; ++iLoc ) {
          // Form global row index from local row index
          const Int i = grid_row_ + iLoc*grid_height_;
          // If diagonal entry, set to one, otherwise zero
          cos_values[iLoc + jLoc * local_height_] = sin(dec_[i] * kDegreeToRadian) *
              sin(dec_[j] * kDegreeToRadian) + cos(dec_[i] * kDegreeToRadian) *
              cos(dec_[j] * kDegreeToRadian) * cos((ra_[i] - ra_[j]) * kDegreeToRadian);
      }
  }

    if (grid_->Rank() == 0) {
      PrintRawArray(cos_values, 10, 1);
    }

  std::vector<double> previous_previous(local_height_ * local_width_, 1.0f);
  std::vector<double> previous = cos_values;
  std::vector<double> current;
  std::vector<double> local_sum = cos_values;
  std::vector<double> local_signal;

  int k = 0;
  local_signal = std::vector<double>(local_height_ * local_width_, 0.0f);
  for (int ell = 1; ell < c_end_[bands_-1]; ++ell) {
    double coefficient = ((double) 2 * ell + 1)/((double) 2 * ell * (ell + 1));
    if (ell == 1){
      current = cos_values;
      VectorTimesScalar(current, coefficient);
      if (ell >= c_start_[k]) VectorPlusEqualsVector(local_signal, current);
      continue;
    }


    current = std::vector<double>(local_height_ * local_width_, 2 * ell - 1);
    VectorTimesEqualsVector(current, cos_values);
    VectorTimesEqualsVector(current, previous);

    VectorTimesScalar(previous_previous, (float)(1 - ell) / (float)ell);

    VectorPlusEqualsVector(current, previous_previous);

    previous_previous = previous;
    previous = current;

    if (ell < c_start_[k]) continue;

    // if (grid_->Rank() == 0) {
    //   std::cout << "Before: Mode " << ell << " coef: " << coefficient << std::endl;
    //   PrintRawArray(previous, 10, 1);
    // }
    VectorTimesScalar(current, coefficient);
    // if (grid_->Rank() == 0) {
    //   std::cout << "After: Mode " << ell << " coef: " << coefficient << std::endl;
    //   PrintRawArray(previous, 10, 1);
    // }

    VectorPlusEqualsVector(local_signal, current);



    
    if (ell>=c_end_[k]){
      signal_[k].Attach(bins_, bins_, *grid_, 0, 0, local_signal.data(), local_height_ );
      //std::cout << "value " << signal_[k].Get(0,0) << std::endl;
      //Print(signal_[k], "Signal");
      // std::cout << "value " << signal_[k].Get(0,0) << std::endl;
      //if (grid_->Rank() == 1) PrintRawArray(current);
      // if (i_have_it) std::cout << "local: " << local_signal[0] << std::endl;
      return;

      ell = c_start_[++k];

      local_signal = std::vector<double>(local_height_ * local_width_, 0.0f);
    }


  }

}

void AngularPowerSpectrum::PrintRawArray(std::vector<double> v, int length, int height) {
  for (int i = 0; i < height; ++i){
    for (int j = 0; j < length; ++j){
      std::cout << " " << v[i + j*height];
    }
    std::cout << std::endl;
  }
}

void AngularPowerSpectrum::SaveDistributedMatrix(std::string name, DistMatrix<double> *matrix, Int num_rows, Int num_cols) {
  std::ofstream outfile (test_directory_ + name, std::ios::binary);
  double number;
  for (int i = 0; i < num_rows; ++i){
    for (int j = 0; j < num_cols; ++j){
      number = matrix->Get(i,j);
      outfile.write(reinterpret_cast<char *>(&number), sizeof(number));
    }
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