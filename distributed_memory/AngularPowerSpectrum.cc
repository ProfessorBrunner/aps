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
       local_width_(Length(bins, grid_col_, grid_width_)),
       is_root_(grid.Rank() == 0) {
}

AngularPowerSpectrum::~AngularPowerSpectrum() {}


void AngularPowerSpectrum::run() {
  Timer timer;
  Barrier();
  timer.Start();
  CalculateSignal();
  Barrier();
  double elapsed = timer.Stop();
  if (is_root_) std::cout << "Signal calculated in " << elapsed << std::endl;
  Print(signal_[0], "Signal_0");
}

/**
 * Build Covariance Matrix
 */
void AngularPowerSpectrum::CalculateCovariance() {}


void AngularPowerSpectrum::CalculateSignal() {
  // TODO(Alex): This code doesn't take advantage of symmetry. It could
  // also calculate each band in i-j blocks to take better advantage of memory
  signal_ = std::vector<DistMatrix<double>>(bands_, DistMatrix<double>(*grid_));
  sum_ = DistMatrix<double>(*grid_);
  //Build Cosine Matrix
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



  //Initialize vectors for Legendre Calculation & DistMatrix operations
  std::vector<double> previous_previous(local_height_ * local_width_, 1.0f);
  std::vector<double> previous = cos_values;
  std::vector<double> current; 
  local_signal = std::vector<std::vector<double>>(bands_, std::vector<double> (local_height_ * local_width_, 0.0f));
  local_sum = std::vector<double>(local_height_ * local_width_, 0.0f);

  //Begin Legendre Calculation
  int k = 0;
  for (int ell = 1; ell <= c_end_[bands_-1]; ++ell) {
    double coefficient = ((double) 2 * ell + 1)/((double) 2 * ell * (ell + 1));
    //Base case: ell = 1
    if (ell == 1){
      current = cos_values;
      VectorTimesScalar(current, coefficient);
      if (ell >= c_start_[k]) VectorPlusEqualsVector(local_signal[k], current);
      continue;
    }
    //current = current * cos_values * previous + (previous_previous * (1-ell/ell))
    current = std::vector<double>(local_height_ * local_width_, ((double) 2 * ell - 1) / (double) ell);
    VectorTimesEqualsVector(current, cos_values);
    VectorTimesEqualsVector(current, previous);
    VectorTimesScalar(previous_previous, (double)(1 - ell) / (double)ell);
    VectorPlusEqualsVector(current, previous_previous);
    //Reset previous_previous and previous for next iteration
    previous_previous = previous;
    previous = current;

    if (ell < c_start_[k]) continue;
    //Set local_signal to current
    VectorTimesScalar(current, coefficient);
    VectorPlusEqualsVector(local_signal[k], current);


    if (ell == c_end_[k]){
      if (is_root_) {
        std::cout << "Attaching band " << k << " modes: " << c_start_[k] <<
            " to " << c_end_[k] << std::endl;
      }
      signal_[k].Attach(bins_, bins_, *grid_, 0, 0, local_signal[k].data(), 
          local_height_ );



      VectorTimesScalar(local_signal[k], c_[k]);
      VectorPlusEqualsVector(local_sum, local_signal[k]);

#     ifdef APS_OUTPUT_TEST
        SaveDistributedMatrix("signal", &signal_[k], bins_, bins_);
#     endif

      ++k;
      //local_signal[k] = std::vector<double>(local_height_ * local_width_, 0.0f);
    }


  }
  sum_.Attach(bins_, bins_, *grid_, 0, 0, local_sum.data(), local_height_ );
  //Print(*sum_, "Sum");
  
}

void AngularPowerSpectrum::PrintRawArray(std::vector<double> v, int length, 
    int height) {
  for (int i = 0; i < height; ++i){
    for (int j = 0; j < length; ++j){
      std::cout << " " << v[i + j*height];
    }
    std::cout << std::endl;
  }
}

void AngularPowerSpectrum::SaveDistributedMatrix(std::string name, 
    DistMatrix<double> *matrix, Int num_rows, Int num_cols) {
  std::ofstream outfile (test_directory_ + name, 
      std::ios::binary | std::ofstream::app );
  double number;
  if (is_root_){
    std::cout << "Writing test file " << name << " to: " << test_directory_ << name << std::endl;
  }
  for (int i = 0; i < num_rows; ++i){
    for (int j = 0; j < num_cols; ++j){
      number = matrix->Get(i,j);
      outfile.write(reinterpret_cast<char *>(&number), sizeof(number));
    }
  }
}


void AngularPowerSpectrum::KLCompression() {

}

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