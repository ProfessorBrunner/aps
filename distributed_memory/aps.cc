/**
 * @file   aps.cc
 * @brief  Calculates Angular Power Spectrum from data files.
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

//#define APS_DEBUG

#include <stdlib.h>

#include <iostream>

#include "elemental-lite.hpp"

#include "OverdensityMap.h"
#include "BandPower.h"
#include "AngularPowerSpectrum.h"

 /*may not need this?
#include ELEM_DIAGONALSCALE_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_INFINITYNORM_INC
#include ELEM_MAXNORM_INC
#include ELEM_ONENORM_INC
#include ELEM_SVD_INC
#include ELEM_UNIFORM_INC*/
using namespace elem;

int main(int argc, char *argv[]) {
  Initialize(  argc, argv  );
  if (argc != 3) {
    std::cout << "Error: Expected 2 arguments" << std::endl
              << "Usage:" << std::endl
              << "./aps [overdensity.fits] [bandpowers.bands]" << std::endl;
    return EXIT_FAILURE;
  }


  int rank = mpi::WorldRank();
  int bands, bins;
  double total_galaxies, omega;

  OverdensityMap *mp;
  BandPower *bp;

  if (mpi::WorldRank() == 0) {
    mp = new OverdensityMap();
    bp = new BandPower();
    mp->LoadFromFile(argv[1]);
    bp->LoadFromFile(argv[2]);
    
    bands = bp->bands_;
    bins = mp->bins_;
    total_galaxies = mp->total_galaxies_;
    omega = mp->omega_;
  }

  //Distribute scalars to all the children
  mpi::Broadcast(&bands, 1, 0, mpi::COMM_WORLD);
  mpi::Broadcast(&bins, 1, 0, mpi::COMM_WORLD);
  mpi::Broadcast(&total_galaxies, 1, 0, mpi::COMM_WORLD);
  mpi::Broadcast(&omega, 1, 0, mpi::COMM_WORLD);


  //set algorithmic blocksize (default 128): SetBlocksize( int blocksize );
  Grid grid( mpi::COMM_WORLD );

  AngularPowerSpectrum aps(bins, bands, total_galaxies, omega, grid);

  if (mpi::WorldRank() == 0) {
    memcpy(aps.c_, bp->c_, bands * sizeof *bp->c_);
    memcpy(aps.c_start_, bp->c_start_, bands * sizeof *bp->c_start_);
    memcpy(aps.c_end_, bp->c_end_, bands * sizeof *bp->c_end_);
    memcpy(aps.ra_, mp->ra_, bins * sizeof *mp->ra_);
    memcpy(aps.dec_, mp->dec_, bins * sizeof *mp->dec_);
  }

  //Distribute arrays to all the children
  mpi::Broadcast(aps.c_, bands, 0, mpi::COMM_WORLD);
  mpi::Broadcast(aps.c_start_, bands, 0, mpi::COMM_WORLD);
  mpi::Broadcast(aps.c_end_, bands, 0, mpi::COMM_WORLD);
  mpi::Broadcast(aps.ra_, bins, 0, mpi::COMM_WORLD);
  mpi::Broadcast(aps.dec_, bins, 0, mpi::COMM_WORLD);

  delete mp;
  delete bp;

  aps.run();

  return EXIT_SUCCESS;
}