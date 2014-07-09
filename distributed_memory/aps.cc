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
#ifdef APS_OUTPUT_TEST
#include <sys/stat.h>
#endif
#include <iostream>

#include "elemental-lite.hpp"

#include "OverdensityMap.h"
#include "BandPower.h"
#include "AngularPowerSpectrum.h"

using namespace elem;

int main(int argc, char *argv[]) {
  Initialize(  argc, argv  );
  if (argc != 3) {
    std::cout << "Error: Expected 2 arguments" << std::endl
              << "Usage:" << std::endl
              << "./aps [overdensity.fits] [bandpowers.bands]" << std::endl;
    return EXIT_FAILURE;
  }

  int bands, bins;
  double total_galaxies, omega;
  std::string output_directory, test_name, test_directory;

  OverdensityMap *mp;
  BandPower *bp;

  //set algorithmic blocksize (default 128): SetBlocksize( int blocksize );
  Grid grid( mpi::COMM_WORLD );

  if (grid.Rank() == 0) {
    mp = new OverdensityMap();
    bp = new BandPower();
    mp->LoadFromFile(argv[1]);
    bp->LoadFromFile(argv[2]);

    bands = bp->bands_;
    bins = mp->bins_;
    total_galaxies = mp->total_galaxies_;
    omega = mp->omega_;

    //Determine output file directories
    output_directory =std::string(argv[2]);
    int char_position = output_directory.find_last_of('/');
    test_name = output_directory.substr(char_position+1);

    output_directory = output_directory.substr(0, char_position);
    
    char_position = test_name.find_last_of('.');
    test_name = test_name.substr(0, char_position);
    test_directory = output_directory + "/" + std::string("test_distributed_") +
        test_name + "/";

    output_directory = output_directory + "/output";


    std::cout << "Data output directory: " << output_directory << std::endl;
#   ifdef APS_OUTPUT_TEST
    std::cout << "Test data output directory: " << test_directory << std::endl;
    mkdir(test_directory.c_str(), 0766);
#   endif
  }

  //Distribute scalars to all the children
  mpi::Broadcast(&bands, 1, 0, mpi::COMM_WORLD);
  mpi::Broadcast(&bins, 1, 0, mpi::COMM_WORLD);
  mpi::Broadcast(&total_galaxies, 1, 0, mpi::COMM_WORLD);
  mpi::Broadcast(&omega, 1, 0, mpi::COMM_WORLD);

  AngularPowerSpectrum aps(bins, bands, total_galaxies, omega, grid);

  if (grid.Rank() == 0) {
    memcpy(aps.c_, bp->c_, bands * sizeof *bp->c_);
    memcpy(aps.c_start_, bp->c_start_, bands * sizeof *bp->c_start_);
    memcpy(aps.c_end_, bp->c_end_, bands * sizeof *bp->c_end_);
    memcpy(aps.ra_, mp->ra_, bins * sizeof *mp->ra_);
    memcpy(aps.dec_, mp->dec_, bins * sizeof *mp->dec_);

    aps.output_directory_ = output_directory;
    aps.test_directory_ = test_directory;
  }

  //Distribute arrays to all the children
  mpi::Broadcast(aps.c_, bands, 0, mpi::COMM_WORLD);
  mpi::Broadcast(aps.c_start_, bands, 0, mpi::COMM_WORLD);
  mpi::Broadcast(aps.c_end_, bands, 0, mpi::COMM_WORLD);
  mpi::Broadcast(aps.ra_, bins, 0, mpi::COMM_WORLD);
  mpi::Broadcast(aps.dec_, bins, 0, mpi::COMM_WORLD);

  if (grid.Rank() == 0) {
    delete mp;
    delete bp;
  }

  aps.run();

  //mpi::Barrier(mpi::COMM_WORLD);
  Finalize();
  return EXIT_SUCCESS;
}