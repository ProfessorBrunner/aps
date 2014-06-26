/**
 * @file   BandPower.cc
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

#include "BandPower.h"

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <string>



BandPower::BandPower() 
    :  bands_(0),
       c_(NULL),
       c_start_(NULL),
       c_end_(NULL){}

BandPower::~BandPower() {
  if (c_) free(c_);
  if (c_start_) free(c_start_);
  if (c_end_) free(c_end_);
}

void BandPower::LoadFromFile(char *file_path) {
  std::cout << std::string(80, '-') << std::endl;
  std::cout << "Loading bandpowers from: " << file_path << std::endl;
  bands_ = CountLines(file_path);

  c_ = (double *) malloc(bands_ * sizeof(double));
  c_start_ = (int *) malloc(bands_ * sizeof(int));
  c_end_ = (int *) malloc(bands_ * sizeof(int));

  std::ifstream band_file;
  std::string dummy;
  band_file.open(file_path);

  int i = 0;
  /// TODO(Alex): Check for errors
  while (band_file >> dummy >> dummy 
         >> c_start_[i] >> c_end_[i] >> c_[i] 
         >> dummy >> dummy) {
    ++i;
    //std::cout << c_start_[i] << " " << c_end_[i] << " " << c_[i] << std::endl;
  }


  std::cout << i << " bands loaded from file with "<< bands_ << " lines" << std::endl;
  band_file.close();
  std::cout << std::string(80, '-') << std::endl;
}



int BandPower::CountLines(char *file_path) {
  int num_lines = 0;
  std::ifstream in(file_path);
  std::string unused;
  
  while ( std::getline(in, unused) ) ++num_lines;

  in.close();
  return num_lines;
}