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

#ifndef APS_OVERDENSITYMAP_H_
#define APS_OVERDENSITYMAP_H_

/**
 *  Overdensity Map. Loads the overdensity map and header data using fitsio and healpix.
 */
class OverdensityMap {
 public:
  float *healpix_map_;     /**< Overdensity map. Density at each pixel. */
  double *ra_;             /**< Right Ascension. Pixel position in astronomical coordinates. */
  double *dec_;            /**< Declination. Pixel position in astronomical coordinates. */
  int nside_;              /**< Subdivisions per side. Healpix level of subdivision. */
  int bins_;               /**< Number of Pixels. Number of pixels in the map. */
  double total_galaxies_;  /**< Number of Galaxies. Total galaxies used to calculate overdensity. */
  double omega_;           /**< Total Area. Total area of usable pixels. */

  /**
   *  A constructor.
   */
  OverdensityMap();
  /**
   *  A destructor.
   */
  ~OverdensityMap();

   /**
    * Load map from a file, sets public member variables.
    */
  void LoadFromFile(
    char *file /**< [in] .fits file to load. */
    );

  /**
   * Load *** from fits header data.
   * 
   * This was based on read_healpix_map in healpix.
   */
  void LoadFitsKeys(
    char *file /**< [in] .fits file to load. */
    );
  
 private:

};

#endif