// Functions that are not used in the distributed
// APS code.


// From angular_power_spectrum.c
int count_usable_pixels(FILE *, int);
long count_cross_Healpix_pixels(char *, float *, float *);
int read_SDSSpix_file(FILE *, double *, double *, double *, double *);
int read_cross_Healpix_file(double *, float *, float *, double *, double *, long);
double calculate_matrices(double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *);
int calculate_overdensity(double *, double *, double *);
int calculate_covariance(double *, double *, double *, double *, double *, double *, double *);
int calculate_cross_covariance(double *, double *, double *, double *, double *, double *, double *);
int calculate_cos_angle(double *, double *, double *);
int read_signal(FILE *, double *);
int calculate_C(double *, double *, double *, double *);

// From math.c
int multiply_diagonal_matrices(double *, double *, double *, int );
double matrix_trace(double *, int );
double multiply_square_matrices(double *, double *, double *, long);
double dot_product(double *, double *, long);
