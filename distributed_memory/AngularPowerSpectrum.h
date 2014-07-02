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

#ifndef APS_ANGULARPOWERSPECTRUM_H_
#define APS_ANGULARPOWERSPECTRUM_H_


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


  AngularPowerSpectrum();
  ~AngularPowerSpectrum();
  
 private:

};

#endif