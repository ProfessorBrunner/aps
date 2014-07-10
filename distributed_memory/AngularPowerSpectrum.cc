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
#include ELEM_HEMM_INC
#include ELEM_ONES_INC
#include ELEM_SYMM_INC
#include ELEM_COPY_INC
//#include ELEM_FUNCS_INC

#include "chealpix.h"
#include "fitsio.h"

#include "AngularPowerSpectrum.h"


using namespace elem;

AngularPowerSpectrum::AngularPowerSpectrum(int bins, int bands, 
    double total_galaxies, double omega, Grid &grid)
    :  bins_(bins),
       bands_(bands),
       total_galaxies_(total_galaxies),
       omega_(omega),
       inverse_density_(omega_/total_galaxies_),
       grid_(&grid),
       c_(new double[bands]),
       c_start_(new int[bands]),
       c_end_(new int[bands]),
       ra_(new double[bins]),
       dec_(new double[bins]),
       local_overdensity_(new double[bins]),
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
  double elapsed;

  CreateOverdensity();

  Barrier();
  timer.Start();
  CalculateSignal();
  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "Signal calculated in " << elapsed << std::endl;
# ifdef APS_OUTPUT_TEST
  for (int k = 0; k < bands_; ++k){
    //Save matrix to file
    SaveDistributedMatrix("signal" + std::to_string(k), &signal_[k]);
  }
# endif

  timer.Start();
  KLCompression();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "KL compression in " << elapsed << std::endl;
}

void AngularPowerSpectrum::CalculateSignal() {
  //Initialize Signal and Sum matrices
  signal_ = std::vector<DistMatrix<double>>(bands_, DistMatrix<double>(*grid_));
  sum_ = DistMatrix<double>(*grid_);

  //Build Cosine Matrix
  /* TODO(Alex): This code doesn't take advantage of symmetry. It could
   * also calculate each band in i-j blocks to take better advantage of memory */
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
      //Attach local_signal data to the distributed signal matrix
      if (is_root_) {
        std::cout << "Attaching band " << k << " modes: " << c_start_[k] <<
            " to " << c_end_[k] << std::endl;
      }
      signal_[k].Attach(bins_, bins_, *grid_, 0, 0, local_signal[k].data(), 
          local_height_ );

      //Sum matrix calculation: current used as a temporary vector with local_signal's data
      current = local_signal[k];
      VectorTimesScalar(current, c_[k]);
      VectorPlusEqualsVector(local_sum, current);

      ++k;
    }
  }
  //Attach local_sum data to the distributed sum matrix
  sum_.Attach(bins_, bins_, *grid_, 0, 0, local_sum.data(), local_height_ );  
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
    DistMatrix<double> *matrix) {
  Write(*matrix, test_directory_ + name, BINARY_FLAT);
}

void AngularPowerSpectrum::KLCompression() {
  DistMatrix<double> temp, B;
  DistMatrix<double> P(*grid_);
  DistMatrix<double,VR,STAR> w(*grid_);
  double p, q;
  int cutoff;
  double *buffer;
  Copy(sum_, temp);
  Ones( P, bins_, bins_);
  //Symm(RIGHT)    C:= a B A    + b C
  //Symm(RIGHT) temp:= p P sum_ + q temp
  p = - kLargeNumber / ( inverse_density_ * (bins_ + inverse_density_) );
  q = 1 / inverse_density_;
  Symm(RIGHT, UPPER, p, sum_, P, q, temp);

  HermitianEig(UPPER, temp, w, B, DESCENDING);
  buffer = B.Buffer();
  for( Int jLoc=0; jLoc<local_width_; ++jLoc ) {
    // Form global column index from local column index
    const Int j = grid_col_ + jLoc*grid_width_;
    for( Int iLoc=0; iLoc<local_height_; ++iLoc ) {
        // Form global row index from local row index
        const Int i = grid_row_ + iLoc*grid_height_;     
        buffer[iLoc+jLoc*local_height_] /= NoiseSqrtAt(i, j);

    }
  }
  MatrixInfo(B);
  for (cutoff = w.Height()-1; cutoff > 0 && w.Get(cutoff, 0) < 1; --cutoff);
  
  View(temp, B, 0, 0, B.Height(), cutoff);

}

void AngularPowerSpectrum::MatrixInfo(DistMatrix<double> &m){
  std::cout << "Width: " << m.Width() << " Height: " << m.Height() << std::endl;
  std::cout << "Total Elements: " << m.Height() * m.Width() << std::endl;
  std::cout << "Memory Elements: " << m.AllocatedMemory() << std::endl;
  std::cout << "Is a view? " << m.Viewing() << std::endl;
}

void AngularPowerSpectrum::MatrixInfo(DistMatrix<double,VR,STAR> &m){
  std::cout << "Width: " << m.Width() << " Height: " << m.Height() << std::endl;
  std::cout << "Total Elements: " << m.Height() * m.Width() << std::endl;
  std::cout << "Memory Elements: " << m.AllocatedMemory() << std::endl;
  std::cout << "Is a view? " << m.Viewing() << std::endl;
}

void AngularPowerSpectrum::CreateOverdensity() {
  overdensity_ = DistMatrix<double, STAR, STAR>(*grid_);
  overdensity_.Attach(bins_, bins_, *grid_, 0, 0, local_overdensity_, bins_);
}

/*********************
 * To Do:
 ********************/
void AngularPowerSpectrum::CalculateCovariance() {}


void AngularPowerSpectrum::EstimateC() {}

void AngularPowerSpectrum::ExpectedCovariance() {}

void AngularPowerSpectrum::CalculateFisher() {}

void AngularPowerSpectrum::RecalculateC_L() {}