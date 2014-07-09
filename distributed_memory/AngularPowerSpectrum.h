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
  /// c_[i] combined coefficient for the band from c_start_[i] to c_end_[i]
  double *c_;
  /// c_start_[i] is the begining of i-th band
  int *c_start_;
  /// c_start_[i] is the end of i-th band
  int *c_end_;
  ///Overdensity vector. Overdensity of pixels at ra_, dec_ 
  double *overdensity_;
  ///Right Ascension. Pixel position in astronomical coordinates.
  double *ra_;
  ///Declination. Pixel position in astronomical coordinates.
  double *dec_;
  ///
  bool is_root_;
  ///
  std::vector<DistMatrix<double>> signal_;
  ///
  DistMatrix<double> sum_;
  ///
  Grid *grid_;
  /// Grid 
  Int grid_height_;
  /// Grid 
  Int grid_width_;
  /// Grid 
  Int grid_row_;
  /// Grid 
  Int grid_col_;
  ///  
  Int local_height_;
  ///  
  Int local_width_;
  ///
  std::vector<std::vector<double>> local_signal;
  ///
  std::vector<double> local_sum;

  ///
  std::string output_directory_;
  ///
  std::string test_directory_;

  ///  Conversion factor for degrees to radians
  static constexpr double kDegreeToRadian = M_PI/180.0;


  ///default constructor
  AngularPowerSpectrum();
  ///non root constructor
  AngularPowerSpectrum(int bins, int bands, double total_galaxies, 
    double omega, Grid &grid);
  ///destructor
  ~AngularPowerSpectrum();
  /**
   * Called from aps' main()
   */
  void run();
  
 private:
  void CalculateCovariance();
  /**
   * Build Signal Matrix
   */

  void CalculateSignal();
  void KLCompression();
  void EstimateC();
  void ExpectedCovariance();
  void CalculateFisher();
  void RecalculateC_L();
  void SaveRawArray();
  void PrintRawArray(std::vector<double> v, int length, int height);
  void SaveDistributedMatrix(std::string name, DistMatrix<double> *matrix, Int num_rows, Int num_cols);

  inline void Barrier(){
    mpi::Barrier(grid_->Comm());
  }

  inline void VectorPlusEqualsVector(std::vector<double> &a, std::vector<double> &b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<double>());
  }

  inline void VectorTimesEqualsVector(std::vector<double> &a, std::vector<double> &b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::multiplies<double>());
  }

  inline void VectorTimesScalar(std::vector<double> &v, double a) {
    std::transform(v.begin(), v.end(), v.begin(), 
        std::bind1st(std::multiplies<double>(), a));
  }


};

#endif