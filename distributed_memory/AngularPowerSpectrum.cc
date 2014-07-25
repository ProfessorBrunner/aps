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

  //Set Up Timers
  double elapsed;
  Timer timer;
  Timer total_timer;
  Timer iteration_timer;
  Barrier();
  total_timer.Start();
  
  CreateOverdensity();

  /*CALCULATE SIGNAL*/
  Barrier();
  timer.Start();

  CalculateSignal();

  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "Signal calculated in " << elapsed << std::endl;

  //Save Signal Matrix to File
# ifdef APS_OUTPUT_TEST
  for (int k = 0; k < bands_; ++k){
    //Save matrix to file
    //TODO(Alex): This could be done with sstream and fails
    //            with more than 999 bands
    std::string file_name("signal");
    if (k<100) file_name += "0";
    if (k<10) file_name += "0";
    file_name += std::to_string(k), signal_[k];
    SaveDistributedMatrix(file_name, signal_[k]);
  }
# endif

  /*KL-COMPRESSION*/
  Barrier();
  timer.Start();

  KLCompression();

  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "KL compression in " << elapsed << std::endl;

# ifdef APS_OUTPUT_TEST
  for (int k = 0; k < bands_; ++k){
    //Save matrix to file
    std::string file_name("kl_signal");
    if (k<100) file_name += "0";
    if (k<10) file_name += "0";
    file_name += std::to_string(k), signal_[k];
    SaveDistributedMatrix(file_name, signal_[k]);
  }
  SaveDistributedMatrix("kl_noise", noise_);
  SaveDistributedMatrix("kl_overdensity", overdensity_);
# endif

  /*CALCULATE DIFFERENCE*/
  Barrier();
  timer.Start();

  CalculateDifference();

  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "Calculate Difference in " << elapsed << std::endl;

  //Save Difference Matrix to File
# ifdef APS_OUTPUT_TEST
    SaveDistributedMatrix("difference" , difference_);
# endif


  /*ITERATIVE ESTIMATION*/
  for (iteration_ = 1; iteration_ <= 3; ++iteration_) {
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


  Barrier();
  total_timer.Stop();
  if (is_root_) std::cout << "Total run time " << elapsed << std::endl;
}





void AngularPowerSpectrum::CreateOverdensity() {
  DistMatrix<double, CIRC, CIRC> temp_ovden(*grid_);
  temp_ovden.Attach(bins_, 1, *grid_, 0, 0, local_overdensity_, bins_);
  overdensity_ = temp_ovden;
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

  /*** Begin Legendre Calculation ***/
  int k = 0;
  for (int ell = 1; ell <= c_end_[bands_-1]; ++ell) {
    double coefficient = ((double) 2 * ell + 1)/((double) 2 * ell * (ell + 1));
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

      //Sum matrix calculation: current used as a temporary vector 
      //with local_signal's data; get's reset at the top of the loop
      current = local_signal[k];
      VectorTimesScalar(current, c_[k]);
      VectorPlusEqualsVector(local_sum, current);

      ++k;
    }
  }
  //Attach local_sum data to the distributed sum matrix
  sum_.Attach(bins_, bins_, *grid_, 0, 0, local_sum.data(), local_height_ );  
}




void AngularPowerSpectrum::KLCompression() {
  is_compressed_ = true;
  DistMatrix<double> temp_eig(*grid_), temp_transform(*grid_), B(*grid_), B_prime(*grid_), P(*grid_);
  DistMatrix<double,VR,STAR> w(*grid_);
  double p, q;
  int cutoff;
  double *buffer;
  Copy(sum_, temp_eig);
  Ones( P, bins_, bins_);

  //initializing p & q to create formula for values of the inverse Noise matrix
  //Noise = kLargeNumber * P + inverse_density_ * I
  //Proof: 
  //math.stackexchange.com/questions/840855/inverse-of-constant-matrix-plus-diagonal-matrix
  p = - kLargeNumber / ( inverse_density_ * (bins_ * kLargeNumber + inverse_density_) );
  q = 1.0 / inverse_density_;

  //Storing Inverse Noise * Sum into temp_eig
  Symm(RIGHT, UPPER, p, sum_, P, q, temp_eig);

  //Read(temp_eig, "distributed_memory/preeigen.dat", BINARY_FLAT);
# ifdef APS_OUTPUT_TEST
  SaveDistributedMatrix("preeigen" , temp_eig);
# endif

  //Elemental Eigensolver: solves temp_eig, eigenvalues > w, eigenvectors > B, 
  //sorted in descending order
  HermitianEig(UPPER, temp_eig, w, B, DESCENDING);

  //Read(B, "distributed_memory/eigenvectors.dat", BINARY_FLAT);
  //Read(w, "distributed_memory/eigenvalues.dat", BINARY_FLAT);
# ifdef APS_OUTPUT_TEST
  SaveDistributedMatrix("eigenvectors" , B);
  DistMatrix<double> temp_w(w);
  SaveDistributedMatrix("eigenvalues" , temp_w);
# endif

//   DistMatrix<double> eigen_diagonal(*grid_);
//   DistMatrix<double> test7(*grid_), test8(*grid_);
//   std::vector<double> diagonal(w.Height());
//   for (int i = 0; i < w.Height(); ++i) diagonal[i] = w.Get(i,0);
//   Diagonal(eigen_diagonal, diagonal);
//   Ones( test7, bins_, bins_);
//   Ones( test8, bins_, bins_);
//   Gemm(NORMAL, NORMAL, 1.0, B, eigen_diagonal, 0.0, test7);
//   Gemm(NORMAL, TRANSPOSE, 1.0, test7, B, 0.0, test8);

// # ifdef APS_OUTPUT_TEST
//   SaveDistributedMatrix("preeigen2" , test8);
// # endif

  //Normalizing: B_{i,j} /= sqrt(N_{ij})
  buffer = B.Buffer();
  for( Int jLoc=0; jLoc<local_width_; ++jLoc ) {
    const Int j = grid_col_ + jLoc*grid_width_;
    for( Int iLoc=0; iLoc<local_height_; ++iLoc ) {
        const Int i = grid_row_ + iLoc*grid_height_;     
        buffer[iLoc+jLoc*local_height_] /= NoiseSqrtAt(i, j);
    }
  }

  //Creating B_prime to keep eigenvectors with eigenvalues > 1
  for (cutoff = w.Height()-1; cutoff > 0 && w.Get(cutoff, 0) < 1; --cutoff);
  bins_ = cutoff + 1;
  View(B_prime, B, 0, 0, B.Height(), bins_);

  //Create new overdensity vector
  DistMatrix<double, VC, STAR>  temp_overdensity_(*grid_);
  Zeros(temp_overdensity_, bins_, 1);
  //Read(temp_overdensity_, "distributed_memory/kl_overdensity.dat", BINARY_FLAT);
  Gemv(TRANSPOSE, 1.0, B_prime, overdensity_, 0.0, temp_overdensity_);
  overdensity_ = temp_overdensity_;

  //Create new Noise and Signal matrices
  Zeros(noise_, bins_, bins_);
  Copy(B, temp_transform);
  //Read(noise_, "distributed_memory/kl_noise.dat", BINARY_FLAT);
  Gemm(TRANSPOSE, NORMAL, kLargeNumber, B_prime, P, inverse_density_, temp_transform);
  Gemm(NORMAL, NORMAL, 1.0, temp_transform, B_prime, 0.0, noise_);
  for (int i = 0; i < bands_; ++i){
    Gemm(TRANSPOSE, NORMAL, 1.0, B_prime, signal_[i], 0.0, temp_transform);
    signal_[i] = DistMatrix<double>(*grid_);
    Zeros( signal_[i], bins_, bins_);
    // TODO(Alex): Check this is correctly resizing local buffer usage
    local_signal[i]  = std::vector<double>();
    Gemm(NORMAL, NORMAL, 1.0, temp_transform, B_prime, 0.0, signal_[i]);
    //Read( signal_[i], "distributed_memory/kl_signal00"+std::to_string(i)+".dat", BINARY_FLAT);
  }
}

void AngularPowerSpectrum::CalculateDifference() {
  if (is_root_) std::cout << "Calculating Difference" << std::endl;
  difference_ = DistMatrix<double>(*grid_);
  local_difference = std::vector<double>(local_height_ * local_width_, 0.0f);
  DistMatrix<double, STAR, STAR> all_overdensity = overdensity_;

  double sum = 0;
  for(int i = 0; i < bins_; ++i) {
    sum += all_overdensity.Get(i, 0);
  }
  std::cout << "Overdensity sum: " << sum << std::endl;

  //Difference = overdensity * overdensity transpose - Noise
  if (!is_compressed_) {
    for( Int jLoc = 0; jLoc < local_width_; ++jLoc ) {
        const Int j = grid_col_ + jLoc*grid_width_;
        for( Int iLoc = 0; iLoc < local_height_; ++iLoc ) {
            const Int i = grid_row_ + iLoc*grid_height_;
            local_difference[iLoc + jLoc * local_height_] = all_overdensity.Get(i, 0) * all_overdensity.Get(j, 0) - NoiseAt(i, j);
        }
    }
  }else{
    double *buffer = noise_.Buffer();
    for( Int jLoc = 0; jLoc < local_width_; ++jLoc ) {
        const Int j = grid_col_ + jLoc*grid_width_;
        for( Int iLoc = 0; iLoc < local_height_; ++iLoc ) {
            const Int i = grid_row_ + iLoc*grid_height_;
            local_difference[iLoc + jLoc * local_height_] = all_overdensity.Get(i, 0) * all_overdensity.Get(j, 0) - buffer[iLoc + jLoc * local_height_];
        }
    }
  }

  difference_.Attach(bins_, bins_, *grid_, 0, 0, local_difference.data(), local_height_ );
}

void AngularPowerSpectrum::EstimateC() {
  DistMatrix<double> P(*grid_), I(*grid_), temp_avg(*grid_);
  std::vector<DistMatrix<double>> A = std::vector<DistMatrix<double>>(bands_, DistMatrix<double>(*grid_));
  Matrix<double> fisher, average, W_prime;

  if (is_root_) std::cout << "Calculating Model Covariance Inverse" << std::endl;

  //Must recalculate sum if compression has occurred
  if (iteration_ != 0 || is_compressed_) {
    if (is_root_) std::cout << "Calculating Sum" << std::endl;
    Zeros( sum_, bins_, bins_);
    for(int k = 0; k < bands_; ++k) {
      Axpy(c_[k], signal_[k], sum_);
    }
  }

  if (!is_compressed_) {
    Ones( P, bins_, bins_);
    Identity( I, bins_, bins_ );
    Axpy(kLargeNumber, P, sum_);
    Axpy(inverse_density_, I, sum_);
  }else{
    Axpy(1.0, noise_, sum_);
  }

  DistMatrix<double>& covariance_inv = sum_;
  SymmetricInverse(LOWER, covariance_inv);

  /* Fisher Matrix Calculation */
  if (is_root_) std::cout << "Calculating Fisher Matrix" << std::endl;
  for (int k = 0; k < bands_; ++k) {
    Zeros(A[k], bins_, bins_);
    /* TODO(Alex): Could use Hemm/Symm for symetric multiplication */
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
  if (is_root_) Print(fisher, "Fisher");


  /* Average Vector Calculation */
  if (is_root_) std::cout << "Calculating Average vector" << std::endl;
  Zeros(temp_avg, bins_, bins_);
  Zeros(average, bands_, 1);
  for (int k = 0; k < bands_; ++k) {
    Gemm(NORMAL, NORMAL, 1.0, A[k], covariance_inv, 0.0, temp_avg);
    average.Set(k, 0, TraceMultiply(difference_, temp_avg));
  }
  
  /* New Window Matrix and C Calculations */
  if (is_root_) {
    std::cout << "Calculating Window Matrix and New C" << std::endl;
    Matrix<double> fisher_inv_sqrt;
    Matrix<double> Y, Y_inv, W, W_prime, Z, temp_Z, row, result;
    std::vector<double> row_sum(bands_, 0.0);
    //Initialize matrices to zero
    Zeros(Y, bands_, bands_);
    Zeros(Z, bands_, bands_);
    Zeros(temp_Z, bands_, bands_);
    Zeros(W_prime, bands_, bands_);

    //Calculate fisher_inv_sqrt
    Copy(fisher, fisher_inv_sqrt);
    SymmetricInverse(LOWER, fisher_inv_sqrt);
    SquareRoot(fisher_inv_sqrt);

    //Calculate Y_inv
    Gemm(NORMAL, NORMAL, 1.0, fisher_inv_sqrt, fisher, 0.0, Y);
    Copy(Y, Y_inv);
    SymmetricInverse(LOWER, Y_inv);

    //Calculate W
    Copy(Y, W);
    for (int j = 0; j < bands_; ++j) {
      for (int i = 0; i < bands_; ++i) {
        row_sum[j] += W.Get(i,j);
      }
      for (int i = 0; i < bands_; ++i) {
        W.Set(i, j, W.Get(i,j) / row_sum[j]);
      }
    }
    
    //Calculate Z & W_prime
    Gemm(NORMAL, NORMAL, 1.0, W, Y_inv, 0.0, temp_Z);
    Gemm(NORMAL, NORMAL, 1.0, temp_Z, fisher_inv_sqrt, 0.0, Z);
    Gemm(NORMAL, NORMAL, 1.0, Z, fisher, 0.0, W_prime);
    
    //Build matrix of C_ls
    Zeros(result, 1, 1);
    for (int i = 0; i < bands_; ++i) {
      View(row, Z, i, 0, 1, bands_);
      Gemm(NORMAL, NORMAL, 1.0, row, average, 0.0, result);
      c_[i] = 0.5 * result.Get(0,0);
    }

    //Save local matrices
# ifdef APS_OUTPUT_TEST
    SaveMatrix("fisher"+std::string("_iter_")+std::to_string(iteration_), fisher);
    SaveMatrix("average"+std::string("_iter_")+std::to_string(iteration_), average);
    SaveMatrix("window"+std::string("_iter_")+std::to_string(iteration_), W_prime);
    Matrix<double> c_matrix;
    c_matrix.Attach(bands_, 1, c_, bands_);
    SaveMatrix("C"+std::string("_iter_")+std::to_string(iteration_), c_matrix);
# endif
  }
  mpi::Broadcast(c_, bands_, 0, grid_->Comm());

  //Save distributed matrix
# ifdef APS_OUTPUT_TEST
  SaveDistributedMatrix("covariance_model"+std::string("_iter_")+std::to_string(iteration_), covariance_inv);
# endif

}


/******************************************************************************
 ******************************* HELPER FUNCTIONS *****************************
 ******************************************************************************/

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
    DistMatrix<double, VC, STAR> &matrix) {
  Write(matrix, test_directory_ + name, BINARY_FLAT);
}

void AngularPowerSpectrum::SaveDistributedMatrix(std::string name, 
    DistMatrix<double> &matrix) {
  Write(matrix, test_directory_ + name, BINARY_FLAT);
}

void AngularPowerSpectrum::SaveMatrix(std::string name, 
    Matrix<double> &matrix) {
  Write(matrix, test_directory_ + name, BINARY_FLAT);
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