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

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <cstring>
#include <string>

#include "chealpix.h"
#include "fitsio.h"
#include "AngularPowerSpectrum.h"

//constructor
AngularPowerSpectrum::AngularPowerSpectrum(){}


//destructor
AngularPowerSpectrum::~AngularPowerSpectrum() {}

//Build Covariance Matrix
void AngularPowerSpectrum::CalculateCovariance() {}

//Build Signal Matrix
void AngularPowerSpectrum::CalculateSignal() {}

//Karhunen-Loeve Compression
void AngularPowerSpectrum::KLCompression() {}