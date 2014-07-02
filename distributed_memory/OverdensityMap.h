/**
 * @file   OverdensityMap.h
 * @brief  Used to load overdensity map and header data from a .fits file.
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

#ifndef APS_OVERDENSITYMAP_H_
#define APS_OVERDENSITYMAP_H_

#define M_PI 3.14159265358979323846

/**
 *  Overdensity Map. Loads the overdensity map and header data using fitsio 
 *  and healpix.
 */
class OverdensityMap {
 public:
  /// Raw healpix map. Overdensity and unused pixels
  float *healpix_map_;
  ///Overdensity vector. Overdensity of pixels at ra_, dec_ 
  double *overdensity_;
  ///Right Ascension. Pixel position in astronomical coordinates.
  double *ra_;
  ///Declination. Pixel position in astronomical coordinates.
  double *dec_;
  ///Subdivisions per side. Healpix level of subdivision.
  int nside_;
  ///Number of Pixels. Number of pixels in the map.
  int bins_;
  ///Number of Galaxies. Total galaxies used to calculate overdensity.
  double total_galaxies_;
  ///Total Area. Total area of usable pixels.
  double omega_;

  ///  Conversion factor for degrees to radians
  static const double kDegreeToRadian = M_PI/180.0;
  ///  Conversion factor for degrees to radians
  static const double kRadianToDegree = 180.0/M_PI;
  /// Square degrees in a sphere
  static const double kSquareDegreePerSphere = 129600.0/M_PI;

  /**
   *  Initialize members to zero.
   */
  OverdensityMap();
  /**
   *  Free all the arrays.
   */
  ~OverdensityMap();

   /**
    * Load map from a file, sets public member variables.
    */
  void LoadFromFile(
    char *file_path /**< [in] .fits file to load. */
    );

  
 private:
  /**
   * Load NSIDE and NGALAXY from fits header data.
   * 
   * This was based on read_healpix_map in healpix.
   */
  void LoadFitsKeys(
    char *file_path /**< [in] .fits file to load. */
    );


  /**
   * Count the number of used pixels and calculate pixel position.
   * 
   * Stores the coordinates in ra_ and dec_ and value in overdensity_.
   * Normalizes the coordinates and stores number of pixels in bins_.
   */
  void ReadHealpixMap();

};

#endif