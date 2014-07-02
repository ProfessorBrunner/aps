/**
 * @file   BandPower.h
 * @brief  Loads initial guess for bandpowers with window start and end.
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

#ifndef APS_BANDPOWER_H_
#define APS_BANDPOWER_H_


class BandPower {
 public:
  /// c_[i] combined coefficient for the band from c_start_[i] to c_end_[i]
  double *c_;
  /// c_start_[i] is the begining of i-th band
  int *c_start_;
  /// c_start_[i] is the end of i-th band
  int *c_end_;
  /// Total number of bands.
  int bands_;

  /**
   *  Initialize members to zero.
   */
  BandPower();

  /**
   *  Free all the arrays.
   */
  ~BandPower();

  /**
    * Load band powers from a file.
    */
  void LoadFromFile(
    char *file_name /**< [in] data file to load. */
    );
  
 private:
  /**
    * Count the number of lines in the file.
    */
  int CountLines(
    char *file_name /**< [in] data file to load. */
    );

};

#endif
