#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include "mkl.h"
#include "omp.h"
#include "chealpix.h"

#define SNR_RETAINED 0.95              // fraction of signal to noise retained.
#define max_iter 3                     // max allowed iterations for C
#define d2r (M_PI/180.0)               // conversion factor for degrees to radians
#define r2d (180.0/M_PI)               // conversion factor for radians to degrees
#define MAXCHARS 10000                 // max allowed characters per object
#define lines_per_object 1             // lines per object
#define contamination 0.0              // stellar contamination
#define num_proc 16                    // number of processors
#define LARGE_NUMBER 1000000.0         // Large number to remove the mean mode in KL compression

// Defined in *_spectrum.c
extern long bins;
extern long bands;
extern double omega;
extern double total_galaxies;

// In tools.c
int object_count(FILE *, int);

// In math.c
double Legendre(double, int);
double trace_multiply(double *, double *, int, int);
int multiply_diagonal_matrices(double *, double *, double *, int );
double matrix_trace(double *, int );
double invert_matrix(double *, long);
double multiply_square_matrices(double *, double *, double *, long);
double multiply_matrices(double *, double *, double *, long, long, long);
double dot_product(double *, double *, long);
void matrix_transpose(double *, double *, long, long);
void matrix_square_root(double *, double *, long);

// In angular_power_spectrum.c
int count_usable_pixels(FILE *, int);
long count_Healpix_pixels(char *, float *);
long count_cross_Healpix_pixels(char *, float *, float *);
int read_bandpower_file(FILE *, double *, int *, int *);
int read_SDSSpix_file(FILE *, double *, double *, double *, double *);
int read_Healpix_file(double *, float *, double *, double *, long);
int read_cross_Healpix_file(double *, float *, float *, double *, double *, long);
int KL_compression(double *, double *, double *, double *, double *, FILE *);
double calculate_matrices(double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *);
double estimate_C(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, double *, int, FILE *, FILE *, FILE *);
int calculate_overdensity(double *, double *, double *);
int calculate_covariance(double *, double *, double *, double *, double *, double *, double *);
int calculate_Healpix_covariance(double *, double *, double *, double *, double *, double *);
int calculate_cross_covariance(double *, double *, double *, double *, double *, double *, double *);
int calculate_cos_angle(double *, double *, double *);
double calculate_signal(double *, double *, int *, int *);
int print_signal(FILE *, double *);
int read_signal(FILE *, double *);
int calculate_difference(double *, double *, double *, double *, double *, double *);
double calculate_products(double *, double *, double *, double *, double *, double *);
int calculate_Fisher(double *, double *);
int print_Fisher(double *, int, FILE *);
int calculate_C(double *, double *, double *, double *);
int calculate_KL_C(double *, double *, double *, double *, FILE *);
int print_values(double *, int *, int *, double *, double *, int, FILE *);
