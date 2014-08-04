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
#include <malloc.h>

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
//#include ELEM_SQUAREROOT_INC //Currently a bug
#include ELEM_AXPY_INC
#include ELEM_HEMM_INC
#include ELEM_ONES_INC
#include ELEM_DIAGONAL_INC
#include ELEM_SYMM_INC
#include ELEM_COPY_INC
#include ELEM_GEMV_INC
#include ELEM_SYMMETRICINVERSE_INC
#include ELEM_INVERSE_INC
#include ELEM_TRACE_INC
#include ELEM_SCALE_INC

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
  if (is_root_){
    std::cout << "Running AngularPowerSpectrum" << std::endl;
    std::cout << "omega/total_galaxies : " << inverse_density_ << std::endl;
    std::cout << "Nodes                : " << grid_->Size() << std::endl;
    std::cout << "Grid height          : " << grid_height_ << std::endl;
    std::cout << "Grid Width           : " << grid_width_ << std::endl;
  }

  //Set Up Timers
  double elapsed;
  Timer timer;
  Timer total_timer;
  Timer iteration_timer;
  Barrier();
  total_timer.Start();
  
  CreateOverdensity();

  // if (is_root_) {
  //   Matrix<double> c_matrix;
  //   c_matrix.Attach(bands_, 1, c_, bands_);
  //   Print(c_matrix, "Original C");
  // }

  /*CALCULATE SIGNAL*/
  Barrier();
  timer.Start();

  CalculateSignal();

  Barrier();
  elapsed = timer.Stop();
  PrintMemory("CalculateSignal-end");
  if (is_root_) std::cout << "TIME: CalculateSignal[] " << elapsed << std::endl;

  //Save Signal Matrix to File
# ifdef APS_OUTPUT_TEST
  for (int k = 0; k < bands_; ++k){
    //Save matrix to file
    //TODO(Alex): This could be done with sstream and fails
    //            with more than 999 bands (not a big deal now)
    std::string file_name("signal");
    if (k<100) file_name += "0";
    if (k<10) file_name += "0";
    file_name += std::to_string(k), signal_[k];
    SaveDistributedMatrix(file_name, signal_[k]);
  }
# endif


# ifdef APS_KL_COMPRESSION
  /*KL-COMPRESSION*/
  PrintMemory("KLCompression-before");
  Barrier();
  timer.Start();

  KLCompression();

  Barrier();
  elapsed = timer.Stop();
  PrintMemory("KLCompression-after");
  if (is_root_) std::cout << "TIME: KLCompression[] " << elapsed << std::endl;

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
# endif

# ifndef APS_KL_COMPRESSION
  /*CREATE NOISE*/
  Barrier();
  timer.Start();
  CreateNoise();
  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "TIME: CreateNoise[] " << elapsed << std::endl;
# endif

  /*CALCULATE DIFFERENCE*/
  Barrier();
  timer.Start();

  CalculateDifference();

  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "TIME: CalculateDifference[] " << elapsed << std::endl;

  //Save Difference Matrix to File
# ifdef APS_OUTPUT_TEST
    SaveDistributedMatrix("difference" , difference_);
# endif


  /*ITERATIVE ESTIMATION*/
  timer.Start();
  if (is_root_) std::cout << "Iterative Estimation";
  for (iteration_ = 1; iteration_ <= 3; ++iteration_) {
    if (is_root_) std::cout << "Iteration " << iteration_ << std::endl;
    PrintMemory("Iteration-" + std::to_string(iteration_));
    Barrier();
    iteration_timer.Start();

    EstimateC();

    Barrier();
    elapsed = iteration_timer.Stop();
    if (is_root_) std::cout << "TIME: EstimateC[iter=" << iteration_ << "] " 
        << elapsed << std::endl;
  }

  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "TIME: TotalIterations[] " << elapsed << std::endl;
  elapsed = total_timer.Stop();
  if (is_root_) std::cout << "TIME: Total[] " << elapsed << std::endl;
}





void AngularPowerSpectrum::CreateOverdensity() {
  if (is_root_) std::cout << "Creating Overdensity" << std::endl;
  DistMatrix<double, CIRC, CIRC> temp_overdensity(*grid_);
  temp_overdensity.Attach(bins_, 1, *grid_, 0, 0, local_overdensity_, bins_);
  overdensity_ = temp_overdensity;
}

void AngularPowerSpectrum::CalculateSignal() {
  if (is_root_) std::cout << "Calculating Signal Matrix" << std::endl;
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
          cos_values[iLoc + jLoc * local_height_] = 
              sin(dec_[i] * kDegreeToRadian) * sin(dec_[j] * kDegreeToRadian)
              + cos(dec_[i] * kDegreeToRadian) *  cos(dec_[j] * kDegreeToRadian)
              * cos((ra_[i] - ra_[j]) * kDegreeToRadian);
      }
  }

  //Initialize vectors for Legendre Calculation & DistMatrix operations
  std::vector<double> previous_previous(local_height_ * local_width_, 1.0f);
  std::vector<double> previous = cos_values;
  std::vector<double> current; 
  local_signal_ = std::vector<std::vector<double>>(bands_, std::vector<double> (local_height_ * local_width_, 0.0f));
  local_sum_ = std::vector<double>(local_height_ * local_width_, 0.0f);

  /*** Begin Legendre Calculation ***/
  int k = 0;
  for (int ell = 1; ell <= c_end_[bands_-1]; ++ell) {
    double coefficient = ((double) 2 * ell + 1)/((double) 2 * ell * (ell + 1));
    if (ell == 1){
      current = cos_values;
      VectorTimesScalar(current, coefficient);
      if (ell >= c_start_[k]) VectorPlusEqualsVector(local_signal_[k], current);
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

    //Set local_signal_ to current
    VectorTimesScalar(current, coefficient);
    VectorPlusEqualsVector(local_signal_[k], current);

    if (ell == c_end_[k]){
      //Attach local_signal_ data to the distributed signal matrix
      if (is_root_) {
        std::cout << "Attaching band " << k << " modes: " << c_start_[k] <<
            " to " << c_end_[k] << std::endl;
      }
      signal_[k].Attach(bins_, bins_, *grid_, 0, 0, local_signal_[k].data(), 
          local_height_ );

      //Sum matrix calculation: current used as a temporary vector 
      //with local_signal_'s data; get's reset at the top of the loop
      current = local_signal_[k];
      VectorTimesScalar(current, c_[k]);
      VectorPlusEqualsVector(local_sum_, current);

      ++k;
    }
  }
  //Attach local_sum_ data to the distributed sum matrix
  sum_.Attach(bins_, bins_, *grid_, 0, 0, local_sum_.data(), local_height_ );  
}




void AngularPowerSpectrum::KLCompression() {
  if (is_root_) std::cout << "Calculating KL Compression" << std::endl;
  is_compressed_ = true;
  DistMatrix<double> temp_eig(*grid_), temp_transform(*grid_), B(*grid_), B_prime(*grid_), P(*grid_);
  DistMatrix<double,VR,STAR> w(*grid_);
  double p, q;
  int cutoff;
  double *buffer;
  double elapsed;
  Timer timer;
  
  if (is_root_) std::cout << "Preparing for eigensolver" << std::endl;
  Barrier();
  timer.Start();
  
  Copy(sum_, temp_eig);
  Ones( P, bins_, bins_);

  //initializing p & q to create formula for values of the inverse Noise matrix
  //Noise = kLargeNumber * P + inverse_density_ * I
  //Proof: 
  //math.stackexchange.com/questions/840855/inverse-of-constant-matrix-plus-diagonal-matrix
  p = - kLargeNumber / ( inverse_density_ * (bins_ * kLargeNumber + inverse_density_) );
  q = 1.0 / inverse_density_;
  // std::cout << "Inverse Noise: " << std::endl;
  // std::cout << "p: " << p << std::endl;
  // std::cout << "q: " << q << std::endl;


  //Storing Inverse Noise * Sum into temp_eig
  Symm(RIGHT, UPPER, p, sum_, P, q, temp_eig);

  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "TIME: KLCompression-prepare[] " << elapsed << std::endl;
  if (is_root_) std::cout << "Calculating Eigenvectors" << std::endl;
  timer.Start();

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

  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "TIME: KLCompression-eigensolve[] " << elapsed << std::endl;
  PrintMemory("KLCompression-before_transform");
  if (is_root_) std::cout << "Transforming matrices" << std::endl;
  timer.Start();

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
  if (is_root_) std::cout << "KL Compression RATIO " << bins_ << ":" << cutoff + 1 << std::endl;
  bins_ = cutoff + 1;
  View(B_prime, B, 0, 0, B.Height(), bins_);

  //Create new overdensity vector
  DistMatrix<double, VC, STAR>  temp_overdensity(*grid_);
  Zeros(temp_overdensity, bins_, 1);
  //Read(temp_overdensity, "distributed_memory/kl_overdensity.dat", BINARY_FLAT);
  Gemv(TRANSPOSE, 1.0, B_prime, overdensity_, 0.0, temp_overdensity);
  overdensity_ = temp_overdensity;

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
    local_signal_[i]  = std::vector<double>();
    Gemm(NORMAL, NORMAL, 1.0, temp_transform, B_prime, 0.0, signal_[i]);
    //Read( signal_[i], "distributed_memory/kl_signal00"+std::to_string(i)+".dat", BINARY_FLAT);
  }

  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "TIME: KLCompression-transform[] " << elapsed << std::endl;
  PrintMemory("KLCompression-end");
}

void AngularPowerSpectrum::CreateNoise() {
  Zeros(noise_, bins_, bins_);

  double *buffer = noise_.Buffer();
  for( Int jLoc = 0; jLoc < local_width_; ++jLoc ) {
        const Int j = grid_col_ + jLoc*grid_width_;
        for( Int iLoc = 0; iLoc < local_height_; ++iLoc ) {
            const Int i = grid_row_ + iLoc*grid_height_;
            buffer[iLoc + jLoc * local_height_] = NoiseAt(i, j);
        }
    }
}

void AngularPowerSpectrum::CalculateDifference() {
  if (is_root_) std::cout << "Calculating Difference" << std::endl;
  difference_ = DistMatrix<double>(*grid_);
  local_difference_ = std::vector<double>(local_height_ * local_width_, 0.0f);
  DistMatrix<double, STAR, STAR> all_overdensity = overdensity_;

  double sum = 0;
  for(int i = 0; i < bins_; ++i) {
    sum += all_overdensity.Get(i, 0);
  }
  if (is_root_) std::cout << "Overdensity sum: " << sum << " ~0" << std::endl;

  //Difference = overdensity * overdensity transpose - Noise
  if (!is_compressed_) {
    for( Int jLoc = 0; jLoc < local_width_; ++jLoc ) {
        const Int j = grid_col_ + jLoc*grid_width_;
        for( Int iLoc = 0; iLoc < local_height_; ++iLoc ) {
            const Int i = grid_row_ + iLoc*grid_height_;
            local_difference_[iLoc + jLoc * local_height_] =
                all_overdensity.Get(i, 0) * all_overdensity.Get(j, 0)
                - NoiseAt(i, j);
        }
    }
  }else{
    double *buffer = noise_.Buffer();
    for( Int jLoc = 0; jLoc < local_width_; ++jLoc ) {
        const Int j = grid_col_ + jLoc*grid_width_;
        for( Int iLoc = 0; iLoc < local_height_; ++iLoc ) {
            const Int i = grid_row_ + iLoc*grid_height_;
            local_difference_[iLoc + jLoc * local_height_] = 
                all_overdensity.Get(i, 0) * all_overdensity.Get(j, 0)
                - buffer[iLoc + jLoc * local_height_];
        }
    }
  }

  difference_.Attach(bins_, bins_, *grid_, 0, 0, local_difference_.data(), local_height_ );
}

void AngularPowerSpectrum::EstimateC() {
  if (is_root_) std::cout << "Estimating C" << std::endl;
  DistMatrix<double> P(*grid_), temp_avg(*grid_);
  std::vector<DistMatrix<double>> A = std::vector<DistMatrix<double>>(bands_, DistMatrix<double>(*grid_));
  Matrix<double> fisher, average, W_prime;
  double elapsed;
  Timer timer;

  if (is_root_) std::cout << "Calculating Model Covariance Inverse" << std::endl;
  Barrier();
  timer.Start();

  //Must recalculate sum if compression has occurred
  if (iteration_ != 0 || is_compressed_) {
    if (is_root_) std::cout << "Calculating Sum" << std::endl;
    Zeros( sum_, bins_, bins_);
    for(int k = 0; k < bands_; ++k) {
      Axpy(c_[k], signal_[k], sum_);
    }
  }

  Axpy(1.0, noise_, sum_);


  //use a different name to make clear that the matrix is different
  DistMatrix<double>& covariance_inv = sum_;
  SymmetricInverse(LOWER, covariance_inv);
  // Read(sum_, "data/test_shared_CL_32_model_4/covariance_model_iter_1.dat", BINARY_FLAT);

  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "TIME: EstimateC-covariance[] " << elapsed << std::endl;
  timer.Start();

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

  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "TIME: EstimateC-fisher[] " << elapsed << std::endl;
  timer.Start();


  /* Average Vector Calculation */
  if (is_root_) std::cout << "Calculating Average vector" << std::endl;
  Zeros(temp_avg, bins_, bins_);
  Zeros(average, bands_, 1);
  PrintMemory("EstimateC-average");
  for (int k = 0; k < bands_; ++k) {
    Gemm(NORMAL, NORMAL, 1.0, A[k], covariance_inv, 0.0, temp_avg);
    average.Set(k, 0, TraceMultiply(difference_, temp_avg));
    A[k].Empty();

  }

  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "TIME: EstimateC-average[] " << elapsed << std::endl;
  timer.Start();
  
  /* New Window Matrix and C Calculations */
  if (is_root_) {
    std::cout << "Calculating Window Matrix and New C" << std::endl;
    Matrix<double> fisher_inv_sqrt;
    Matrix<double> fisher_inv;
    Matrix<double> Y, Y_inv, W, W_prime, Z, temp_Z, row, result;
    //Initialize matrices to zero
    Zeros(Y, bands_, bands_);
    Zeros(Z, bands_, bands_);
    Zeros(temp_Z, bands_, bands_);
    Zeros(W_prime, bands_, bands_);
    Zeros(fisher_inv_sqrt, bands_, bands_);


    //Calculate fisher_inv_sqrt
    Copy(fisher, fisher_inv);
    SymmetricInverse(LOWER, fisher_inv);

    // Read(fisher_inv, "data/test_shared_CL_32_model_4/inv_fisher_iter_1.dat", BINARY_FLAT);

    Copy(fisher_inv, fisher_inv_sqrt);
    SquareRoot(fisher_inv_sqrt); //Elemental version has bug


    // Matrix<double> test;
    // Zeros(test, bands_, bands_);
    // Gemm(NORMAL, NORMAL, 1.0, fisher_inv_sqrt, fisher_inv_sqrt, 0.0, test);
    // Gemm(NORMAL, NORMAL, 1.0, test, fisher, 0.0 test);
    // Print(test, "Test Fisher sqrt back to inverse");

    //Calculate Y_inv
    Gemm(NORMAL, NORMAL, 1.0, fisher_inv_sqrt, fisher, 0.0, Y);
    Copy(Y, Y_inv);
    Inverse(Y_inv);

    //Calculate W
    Copy(Y, W);
    for (int i = 0; i < bands_; ++i) {
      double sum = 0.0;
      for (int j = 0; j < bands_; ++j) {
        sum += W.Get(i,j);
      }
      for (int j = 0; j < bands_; ++j) {
        W.Set(i, j, W.Get(i,j) / sum);
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

    //Save bands into output file
    FILE* band_file;
    std::string band_file_name = output_directory_ + std::string("/C_") + 
        output_name_ + std::string(".bands");
    if (iteration_ == 1) {
      band_file = fopen(band_file_name.c_str(), "w");
    }else{
      band_file = fopen(band_file_name.c_str(), "a");
    }
    for (int i = 0; i < bands_; ++i) {
      fprintf(band_file, "%d %d %d %d %lf %lf %lf\n",
          i, (c_end_[i]+c_start_[i])/2, c_start_[i], c_end_[i],
          c_[i], sqrt(fisher_inv.Get(i,i)), c_[i]);
    }
    fclose(band_file);

    //Save local matrices
# ifdef APS_OUTPUT_TEST
    SaveMatrix(std::string("inv_sqrt_fisher_iter_")+std::to_string(iteration_), fisher_inv_sqrt);
    SaveMatrix(std::string("inv_fisher_iter_")+std::to_string(iteration_), fisher_inv);
    SaveMatrix(std::string("fisher_iter_")+std::to_string(iteration_), fisher);
    SaveMatrix(std::string("average_iter_")+std::to_string(iteration_), average);
    SaveMatrix(std::string("window_iter_")+std::to_string(iteration_), W_prime);
    SaveMatrix(std::string("pre_window_iter_")+std::to_string(iteration_), W);
    SaveMatrix(std::string("Y_iter_")+std::to_string(iteration_), Y_inv);
    Matrix<double> c_matrix;
    c_matrix.Attach(bands_, 1, c_, bands_);
    SaveMatrix(std::string("C_iter_")+std::to_string(iteration_), c_matrix);
# endif
  }
  Barrier();
  elapsed = timer.Stop();
  if (is_root_) std::cout << "TIME: EstimateC-c[] " << elapsed << std::endl;
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

void AngularPowerSpectrum::PrintMemory(std::string location) {
  struct mallinfo m;
  m = mallinfo();
  for (int i = 0; i<grid_->Size(); ++i) {
    if (grid_->Rank() == i) std::cout << "MALLOC: " << location
              << "." << i
              << " total: " << m.hblkhd + m.uordblks
              << " mmap: " << m.hblkhd
              << " chunks: " << m.uordblks
              << " sbrk: " << m.arena << std::endl;
    Barrier();
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


void AngularPowerSpectrum::SquareRoot(Matrix<double> &m) {
  Matrix<double> eigen_vectors, eigen_values, temp, sqrt_diagonal;

  Copy(m, temp);
  HermitianEig(UPPER, temp, eigen_values, eigen_vectors, ASCENDING);

  int n = eigen_values.Height();
  std::vector<double> diagonal(n);
  for (int i = 0; i < n; ++i) diagonal[i] = sqrt(eigen_values.Get(i,0));
  Diagonal(sqrt_diagonal, diagonal);

  Gemm(NORMAL, NORMAL, 1.0, eigen_vectors, sqrt_diagonal, 0.0, temp);
  Gemm(NORMAL, TRANSPOSE, 1.0, temp, eigen_vectors, 0.0, m);
}