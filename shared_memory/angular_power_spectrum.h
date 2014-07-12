#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include "mkl.h"
#include "omp.h"
#include "chealpix.h"

#define kSnrRetained 0.95              // fraction of signal to noise retained.
#define kMaxIter 3                     // max allowed iterations for C
#define kDegreeToRadian (M_PI/180.0)   // conversion factor for degrees to radians
#define kRadianToDegree (180.0/M_PI)   // conversion factor for radians to degrees
#define kMaxChars 10000                // max allowed characters per object
#define kLinesPerObject 1              // lines per object
#define kContamination 0.0             // stellar contamination
#define kNumThreads 16                 // number of processors
#define KLargeNumber 1000000.0         // Large number to remove the mean mode in KL compression
//#define APS_OUTPUT_TEST              // enabled through gcc flag
//#define APS_SUPPRESS_MATRIX_STDOUT   // enabled through gcc flag

/* Defined in *_spectrum.c */
// Total number of usable pixels.
extern long g_bins;
// Total number of bandpowers.
extern long g_bands;
// Total area of usable pixels.
extern double g_omega;
// Total number of galaxies in usable pixels.
extern double g_total_galaxies;

/*In tools.c*/
typedef struct Timer_s {
  struct timeval start;
  struct timeval stop;
  double elapsed;
} Timer;

//count number of lines in a file
int object_count(FILE *, int);
void tic(Timer *timer); //Start timer
void toc(Timer *timer); //Stop timer
void toc_print(Timer *timer); //Stop timer and print time
//Write an array to disk in file "dir/name"
void save_raw_double_array(char *root, char *name, double *data, int len);
/* In math.c*/


/* calculate_signal: (cos_angle[j], k)*/
double Legendre(double, int); 
/*
 * calculate_products: (difference, &B[l*bins*bins], bins, bins)
 * calculate_Fisher: (&A[l*bins*bins], &A[l_prime*bins*bins], bins, bins)
 */
double trace_multiply(double *, double *, int, int); 
/* Not used ever */
int multiply_diagonal_matrices(double *, double *, double *, int );
/* Not used ever */
double matrix_trace(double *, int );
/*
 * KL_compression: (noise, bins) x 2
 * estimate_C: (model_covariance, bins), (F, bands) 
 * calculate_KL_C: invert_matrix(A, bands) x 2, (F, bands) x 2
 */
double invert_matrix(double *, long);
/* Not used ever */
double multiply_square_matrices(double *, double *, double *, long);
/*
 * KL_compression: (B_primet, overdensity, test1, k, bins, 1), (B_primet, noise, test5, k, bins, bins), 
 *	(test5, B_prime, test6, k, bins, k), (B_primet, &signal[l*bins*bins], test5, k, bins, bins), 
 *	(test5, B_prime, &test7[l*k*k], k, bins, k), 
 * calculate_KL_C: (invsqrtF, F, A, bands, bands, bands), (W, A, D, bands, bands, bands), 
 * 	(D, invsqrtF, B, bands, bands, bands), (B, F, W, bands, bands, bands)
 */
double multiply_matrices(double *, double *, double *, long, long, long);
/* Not used ever */
double dot_product(double *, double *, long);
/* 
 * KL_compression: (B, Bt, bins, bins), (B_prime, B_primet, bins, k)
 */
void matrix_transpose(double *, double *, long, long);
/* calculate_Fisher: (F, invsqrtF, bands) */
void matrix_square_root(double *, double *, long);


/* 
 * In angular_power_spectrum.c 
 */


long count_Healpix_pixels(char *, float *);
int read_bandpower_file(FILE *, double *, int *, int *);
int read_Healpix_file(double *, float *, double *, double *, long);
int KL_compression(double *, double *, double *, double *, double *, FILE *, char* test_root);
double estimate_C(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, double *, int, FILE *, FILE *, FILE *);
int calculate_Healpix_covariance(double *, double *, double *, double *, double *, double *);
double calculate_signal(double *, double *, int *, int *);
int print_signal(FILE *, double *);
int calculate_difference(double *, double *, double *, double *, double *, double *);
double calculate_products(double *, double *, double *, double *, double *, double *);
int calculate_Fisher(double *, double *);
int print_Fisher(double *, int, FILE *);
int calculate_KL_C(double *, double *, double *, double *, FILE *);
int print_values(double *, int *, int *, double *, double *, int, FILE *);