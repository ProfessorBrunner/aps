/**
 * @file   AngularPowerSpectrum.h
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

#ifndef APS_ANGULARPOWERSPECTRUM_H_
#define APS_ANGULARPOWERSPECTRUM_H_

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include <vector>
#include <stdlib.h>

#include "elemental-lite.hpp"

using namespace elem;

class AngularPowerSpectrum {
 public:
  ///Number of Pixels. Number of pixels in the map.
  int bins_;
  /// Total number of bands.
  int bands_;
  ///Number of Galaxies. Total galaxies used to calculate overdensity.
  double total_galaxies_;
  ///Total Area. Total area of usable pixels.
  double omega_;
  ///omega_/total_galaxies_
  double inverse_density_;
  /// c_[i] combined coefficient for the band from c_start_[i] to c_end_[i]
  double *c_;
  /// c_start_[i] is the begining of i-th band
  int *c_start_;
  /// c_start_[i] is the end of i-th band
  int *c_end_;
  ///Overdensity vector. Overdensity of pixels at ra_, dec_ 
  DistMatrix<double, VC, STAR> overdensity_;
  ///The local copy of Overdensity
  double *local_overdensity_;
  ///Right Ascension. Pixel position in astronomical coordinates.
  double *ra_;
  ///Declination. Pixel position in astronomical coordinates.
  double *dec_;
  /// A boolean check if this process is the root.
  bool is_root_;
  /// A vector of DistMatrices representing the signal matrix.
  std::vector<DistMatrix<double>> signal_;
  /// Sum matrix. Calculated in calculate_signal. Used in KL-compression.
  DistMatrix<double> sum_;
  /// Communication grid. Used for Distributed Matrix functions.
  Grid *grid_;
  /// Grid Height
  Int grid_height_;
  /// Grid Width
  Int grid_width_;
  /// Grid Row Align
  Int grid_row_;
  /// Grid Column Align
  Int grid_col_;
  ///  The local height of a DistMatrix
  Int local_height_;
  ///  The local width of a DistMatrix
  Int local_width_;
  ///  The local slice of the signal matrix in a vector of vectors.
  std::vector<std::vector<double>> local_signal;
  ///  The local slice of the sum matrix.
  std::vector<double> local_sum;
  ///  Output directory for writing
  std::string output_directory_;
  ///  Test directory for writing
  std::string test_directory_;

  ///  Conversion factor for degrees to radians
  static constexpr double kDegreeToRadian = M_PI/180.0;
  ///  Large number to eliminate mean mode in KL compression
  static constexpr double kLargeNumber = 1000000; 


  ///default constructor
  AngularPowerSpectrum();
  ///non root constructor
  AngularPowerSpectrum(int bins, int bands, double total_galaxies, 
    double omega, Grid &grid);
  ///destructor
  ~AngularPowerSpectrum();

  /**
   * Called from aps' main()
   * Root function that calls all other functions.
   */
  void run();
  
 private:

  void MatrixInfo(DistMatrix<double> &m);
  void MatrixInfo(DistMatrix<double,VR,STAR> &m);

  void CreateOverdensity();
  /**
   * Build Covariance Matrix
   */
  void CalculateCovariance();

  /**
   * Build Signal and Sum Matrices
   */
  void CalculateSignal();

  /**
   * Perform Karhunen-Loeve Compression on the data vectors and matrices
   */
  void KLCompression();

  /**
   * Estimates the APS given Signal & Covariance Matrices
   */
  void EstimateC();

  /**
   * Build Expected Covariance Matrix based on C_l
   * Previously named calculate_difference 
   */
  void ExpectedCovariance();

  /**
   * Build Fisher Matrix & Weighted Average
   */
  void CalculateFisher();

  /**
   * Recalculate C_l from the Fisher Matrix
   */
  void RecalculateC_L();

  /**
   * Prints the given vector
   */
  void PrintRawArray(std::vector<double> v, int length, int height);

  /*
   * Saves the DistMatrix with the given file name using Elemental's Write function
   */
  void SaveDistributedMatrix(std::string name, DistMatrix<double> *matrix);

  /// Local barrier method
  inline void Barrier(){
    mpi::Barrier(grid_->Comm());
  }

  /// Local += method for vectors
  inline void VectorPlusEqualsVector(std::vector<double> &a, std::vector<double> &b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<double>());
  }

  ///Local *= method for vectors
  inline void VectorTimesEqualsVector(std::vector<double> &a, std::vector<double> &b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::multiplies<double>());
  }

  ///Local vector  *= scalar method
  inline void VectorTimesScalar(std::vector<double> &v, double a) {
    std::transform(v.begin(), v.end(), v.begin(), 
        std::bind1st(std::multiplies<double>(), a));
  }

  inline double NoiseSqrtAt(int i, int j) {
    return sqrt(((double) i == j) * inverse_density_ + kLargeNumber); 
  }
};
#endif