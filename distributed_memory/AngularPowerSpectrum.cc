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
#include <stdexcept>

#include "elemental-lite.hpp"
#include ELEM_IDENTITY_INC
#include ELEM_DOT_INC
#include ELEM_SQUAREROOT_INC
#include ELEM_AXPY_INC
#include ELEM_HEMM_INC
#include ELEM_ONES_INC
#include ELEM_DIAGONAL_INC
#include ELEM_SYMM_INC
#include ELEM_COPY_INC
#include ELEM_GEMV_INC
#include ELEM_SYMMETRICINVERSE_INC
#include ELEM_TRACE_INC

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
       inverse_density_(omega/total_galaxies),
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
       is_root_(grid.Rank() == 0),
       is_compressed_(false) {
}

AngularPowerSpectrum::~AngularPowerSpectrum() {}





void AngularPowerSpectrum::run() {
  if (is_root_) std::cout << "Running AngularPowerSpectrum" << std::endl;
  Timer timer;
  Timer iteration_timer;
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
    SaveDistributedMatrix("signal" + std::to_string(k), signal_[k]);
  }
# endif

  // Barrier();
  // timer.Start();
  // KLCompression();
  // Barrier();
  // elapsed = timer.Stop();
  // if (is_root_) std::cout << "KL compression in " << elapsed << std::endl;
  // # ifdef APS_OUTPUT_TEST
  // for (int k = 0; k < bands_; ++k){
  //   //Save matrix to file
  //   SaveDistributedMatrix("kl_signal" + std::to_string(k), signal_[k]);
  // }
  //# endif

  Barrier();
  timer.Start();
  //Print(overdensity_, "overdensity_");
  CalculateDifference();
  for (iteration_ = 1; iteration_ <= 1; ++iteration_) {
    Barrier();
    iteration_timer.Start();
    EstimateC();
    Barrier();
    elapsed = iteration_timer.Stop();
    if (is_root_) std::cout << "Iteration " << iteration_ << " time: "<< elapsed << std::endl;
  }
  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "Total iteration time " << elapsed << std::endl;
# ifdef APS_OUTPUT_TEST
    SaveDistributedMatrix("difference" , difference_);
# endif
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
    DistMatrix<double> &matrix) {
  Write(matrix, test_directory_ + name, BINARY_FLAT);
}

void AngularPowerSpectrum::SaveMatrix(std::string name, 
    Matrix<double> &matrix) {
  Write(matrix, test_directory_ + name, BINARY_FLAT);
}

void AngularPowerSpectrum::KLCompression() {
  is_compressed_ = true;
  DistMatrix<double> temp(*grid_), B(*grid_), B_prime(*grid_);
  DistMatrix<double> P(*grid_);
  DistMatrix<double> test7(*grid_);
  DistMatrix<double,VR,STAR> w(*grid_);
  double p, q;
  int cutoff;
  double *buffer;
  Copy(sum_, temp); //Why is this?
  Ones( P, bins_, bins_);
  p = - kLargeNumber / ( inverse_density_ * (bins_ * kLargeNumber + inverse_density_) );
  q = 1.0 / inverse_density_;
  //shared
  //noise: 1000000.765106 1000000
  //inverse noise: 1.30531 -0.00170183

  //distributed
  //p -1700.14 q 1.30701 inverse density 0.765106
  //inverse noise -1698.83 -1700.14
  //std::cout << "p " << p << " q " << q << " inverse density " << inverse_density_<<std::endl;
  //std::cout << "inverse noise "<< p + q << " " << p << std::endl;
  //Symm(RIGHT)    C:= a B A    + b C
  //Symm(RIGHT) temp:= p P sum_ + q temp
  Symm(RIGHT, UPPER, p, sum_, P, q, temp);

  //Print(temp, "Vector to be eigenvalued");

  HermitianEig(UPPER, temp, w, B, DESCENDING);

  // DistMatrix<double> eigen_diagonal(*grid_);
  // DistMatrix<double> test7(*grid_), test8(*grid_);
  // std::vector<double> diagonal(w.Height());
  // for (int i = 0; i < w.Height(); ++i) diagonal[i] = w.Get(i,0);
  // Diagonal(eigen_diagonal, diagonal);
  // Ones( test7, bins_, bins_);
  // Ones( test8, bins_, bins_);
  // Gemm(NORMAL, NORMAL, 1.0, B, eigen_diagonal, 0.0, test7);
  // Gemm(NORMAL, TRANSPOSE, 1.0, test7, B, 0.0, test8);
  // Print(B, "Eigenvectors");
  // Print(eigen_diagonal, "Diagonal");
  // Print(temp, "Original Noise*Sum");
  // Print(test8, "Recreated");
  

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









  SaveDistributedMatrix("eigenvectors", B);
  //MatrixInfo(B);

  for (cutoff = w.Height()-1; cutoff > 0 && w.Get(cutoff, 0) < 1; --cutoff);
  bins_ = cutoff + 1;
  //B_prime is a view of B
  View(B_prime, B, 0, 0, B.Height(), bins_);

  auto temp_overdensity_(overdensity_);
  Gemv(TRANSPOSE, 1.0, B_prime, overdensity_, 0.0, temp_overdensity_);
  overdensity_ = temp_overdensity_;
  Print(overdensity_, "overdensity");
  // Print(w, "eigenvalues");

  //Noise matrix is kLargeNumber*P+inverse_density_*I
  Copy(B, test7);
  Gemm(TRANSPOSE, NORMAL, kLargeNumber, B_prime, P, inverse_density_, test7);
  Gemm(NORMAL, NORMAL, 1.0, test7, B_prime, 0.0, noise_);
  for (int i = 0; i < bands_; ++i){
    Gemm(TRANSPOSE, NORMAL, 1.0, B_prime, signal_[i], 0.0, test7);
    signal_[i] = DistMatrix<double>(*grid_);
    Zeros( signal_[i], bins_, bins_);
    local_signal[i]  = std::vector<double>();
    Gemm(NORMAL, NORMAL, 1.0, test7, B_prime, 0.0, signal_[i]);
  }
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
  DistMatrix<double, CIRC, CIRC> temp(*grid_);
  temp.Attach(bins_, 1, *grid_, 0, 0, local_overdensity_, bins_);
  overdensity_ = temp;
}

void AngularPowerSpectrum::CalculateDifference() {
  difference_ = DistMatrix<double>(*grid_);
  local_difference = std::vector<double>(local_height_ * local_width_, 0.0f);

  for( Int jLoc = 0; jLoc < local_width_; ++jLoc ) {
      // Form global column index from local column index
      const Int j = grid_col_ + jLoc*grid_width_;
      for( Int iLoc = 0; iLoc < local_height_; ++iLoc ) {
          // Form global row index from local row index
          const Int i = grid_row_ + iLoc*grid_height_;
          // If diagonal entry, set to one, otherwise zero
          local_difference[iLoc + jLoc * local_height_] = overdensity_.Get(i, 0) * overdensity_.Get(j, 0) - NoiseAt(i, j);
      }
  }
  difference_.Attach(bins_, bins_, *grid_, 0, 0, local_difference.data(), local_height_ );
}

void AngularPowerSpectrum::EstimateC() {
  DistMatrix<double> P(*grid_), I(*grid_), temp1(*grid_), temp2(*grid_);
  std::vector<DistMatrix<double>> A = std::vector<DistMatrix<double>>(bands_, DistMatrix<double>(*grid_));
  Matrix<double> fisher, average;

  if (iteration_ != 0 || is_compressed_) {
    //Must recalculate sum
    Zeros( sum_, bins_, bins_);
    for(int k = 0; k < bands_; ++k) {
      Axpy(c_[k], signal_[k], sum_);
    }
  }

  if (is_root_) std::cout << "Calculating Model Covariance Inverse" << std::endl;
  Ones( P, bins_, bins_);
  Identity( I, bins_, bins_ );
  Axpy(kLargeNumber, P, sum_);
  Axpy(inverse_density_, I, sum_);

  DistMatrix<double>& covariance_inv = sum_; //more apt variable name
  SymmetricInverse(LOWER, covariance_inv); //this didn't work with UPPER

  /* TODO(Alex): Could use Hemm for symetric multiplication */
  if (is_root_) std::cout << "Fisher Matrix" << std::endl;
  for (int k = 0; k < bands_; ++k) {
    Zeros(A[k], bins_, bins_);
    Gemm(NORMAL, NORMAL, 1.0, covariance_inv, signal_[k], 0.0, A[k]);
  }

  Zeros(fisher, bands_, bands_);
  for (int k = 0; k < bands_; ++k) {
    for (int k_p = k; k_p < bands_; ++k_p) {
      double result = 0.5 * TraceMultiply(A[k], A[k_p]);
      fisher.Set(k_p, k, result);
      fisher.Set(k, k_p, result);
    }
  }
  Print(fisher, "fisher");

  if (is_root_) std::cout << "Calculating Average vector" << std::endl;
  Zeros(temp1, bins_, bins_);
  Zeros(temp2, bins_, bins_);
  Zeros(average, bands_, 1);
  for (int k = 0; k < bands_; ++k) {
    Gemm(NORMAL, NORMAL, 1.0, A[k], covariance_inv, 0.0, temp1);
    average.Set(k, 0, TraceMultiply(difference_, temp1));
  }
  

  if (is_root_) {
    std::cout << "Calculating Window Matrix and New C" << std::endl;
    Matrix<double> fisher_inv_sqrt;
    Matrix<double> Y, Y_inv, W, Z, temp, W_prime, row, result;
    std::vector<double> row_sum(bands_, 0.0);

    Copy(fisher, fisher_inv_sqrt);
    SymmetricInverse(LOWER, fisher_inv_sqrt);
    SquareRoot(fisher_inv_sqrt);

    Zeros(Y, bands_, bands_);

    Gemm(NORMAL, NORMAL, 1.0, fisher_inv_sqrt, fisher, 0.0, Y);
    Copy(Y, Y_inv);
    SymmetricInverse(LOWER, Y_inv);

    Copy(Y, W);
    for (int j = 0; j < bands_; ++j) {
      for (int i = 0; i < bands_; ++i) {
        row_sum[j] += W.Get(i,j);
      }
      for (int i = 0; i < bands_; ++i) {
        W.Set(i, j, W.Get(i,j) / row_sum[j]);
      }
    }
    
    Zeros(Z, bands_, bands_);
    Zeros(temp, bands_, bands_);
    Zeros(W_prime, bands_, bands_);
    Gemm(NORMAL, NORMAL, 1.0, W, Y_inv, 0.0, temp);

    Gemm(NORMAL, NORMAL, 1.0, temp, fisher_inv_sqrt, 0.0, Z);
    Gemm(NORMAL, NORMAL, 1.0, Z, fisher, 0.0, W_prime);
    
    Zeros(result, 1, 1);
    for (int i = 0; i < bands_; ++i) {
      View(row, Z, i, 0, 1, bands_);
      Gemm(NORMAL, NORMAL, 1.0, row, average, 0.0, result);
      c_[i] = 0.5 * result.Get(0,0);
    }
  }
# ifdef APS_OUTPUT_TEST
  SaveDistributedMatrix(std::string("iter_")+std::to_string(iteration_)+"_covariance_model" , covariance_inv);
  if (is_root_) {
    SaveMatrix(std::string("iter_")+std::to_string(iteration_)+"_fisher" , fisher);
    SaveMatrix(std::string("iter_")+std::to_string(iteration_)+"_average" , average);
    Matrix<double> c_matrix;
    c_matrix.Attach(bands_, 1, c_, bands_);
    SaveMatrix(std::string("iter_")+std::to_string(iteration_)+"_C" , c_matrix);
  }
# endif

}

double AngularPowerSpectrum::TraceMultiply(DistMatrix<double> &m1, DistMatrix<double> &m2) {
  DistMatrix<double> row, col, result;
  double sum = 0;
  Int size = m1.Height();

  if (size != m1.Width() || m2.Width() != m2.Height() || size != m2.Width()) {
    //Check if everything is square and same
    throw std::logic_error("Trace Multiply Dimension mismatch");
  }

  Zeros(result, 1, 1);
  for(Int i = 0; i < size; ++i) {
    View(row, m1, i, 0, 1, size);
    View(col, m2, 0, i, size, 1);
    Gemm(NORMAL, NORMAL, 1.0, row, col, 0.0, result);
    //sum += Dot(row, col); //complains about matrix dimensions
    sum += result.Get(0,0);
  }

  return sum;
}