#include "angular_power_spectrum.h"
//#include "new_angular_power_spectrum.h"

/// A function to count the pixels with non-zero area in the objects file.  
int count_usable_pixels(FILE *objects, int number_lines)
{
  int i;
  double lam, et, trash1, trash2, pixel_area;
  char line[MAXCHARS];

  bins = 0;
  for (i=0; i<number_lines; i++) {
    if (fgets(line, sizeof(line), objects)==NULL) printf("Blank line. %s\n", line);
    if (line[0]!='#') {
      if (sscanf(line, "%lf %lf %lf %lf %lf", &trash1, &trash2, &et, &lam, &pixel_area)==5) {
	if (trash2 != 0.0) {
	  bins++;
	}
      } else {printf("Error1: %s\n", line);}
    } else {i--;}
  }
  rewind(objects);
  printf("#Done counting %ld pixels.\n", bins);

  return 0;
}

/// Open the Healpix file and count the pixels and galaxies.
long count_Healpix_pixels(char *input_healpix, float *healpix)
{
  int i;
  long NSIDE;
  char *trash[MAXCHARS], name[MAXCHARS]; // ordering[10], coords[1];

  // Strip the directory structure from the filename if present
  if (strrchr(input_healpix, '/')!=NULL) {
    strcpy(name, strrchr(input_healpix, '/')+1);
  } else {
    strcpy(name, input_healpix);
  }

  NSIDE = strtol(strpbrk(name, "0123456789"), trash, 10);
  total_galaxies = strtod(strpbrk(*trash, "0123456789"), (char**)NULL);
  printf("#Done reading %ld lines.\n", nside2npix(NSIDE));

  bins = 0;
  for (i=0; i<nside2npix(NSIDE); i++) {
    if (healpix[i] >= -1.0) {
      bins++;
    }
  }
  printf("#Done counting %ld pixels.\n", bins);

  return NSIDE;
}

/// Open the Healpix file and count the pixels and galaxies.
long count_cross_Healpix_pixels(char *input_healpix, float *healpix_1, float *healpix_2)
{
  int i;
  long NSIDE;
  char *trash[MAXCHARS], name[MAXCHARS]; // ordering[10], coords[1];

  // Strip the directory structure from the filename if present
  if (strrchr(input_healpix, '/')!=NULL) {
    strcpy(name, strrchr(input_healpix, '/')+1);
  } else {
    strcpy(name, input_healpix);
  }

  NSIDE = strtol(strpbrk(name, "0123456789"), trash, 10);
  total_galaxies = strtod(strpbrk(*trash, "0123456789"), (char**)NULL);
  printf("#Done reading %ld lines.\n", nside2npix(NSIDE));

  bins = 0;
  for (i=0; i<nside2npix(NSIDE); i++) {
    if ((healpix_1[i]>=-1.0) && (healpix_2[i]>=-1.0)) {
      bins++;
    }
  }
  printf("#Done counting %ld pixels.\n", bins);

  return NSIDE;
}

/// A function to read the bandpowers file, and store the results in the appropriate arrays. 
int read_bandpower_file(FILE *bandpowers, double *C, int *C_start, int *C_stop)
{
  int i, trash, garbage;
  double trash1, trash2;
  char line[MAXCHARS];

  for (i=0; i<bands; i++) {
    if (fgets(line, sizeof(line), bandpowers)==NULL) printf("Blank line. %s\n", line);
    if (line[0]!='#') {
      if (sscanf(line, "%d %d %d %d %lf %lf %lf", &trash, &garbage, &C_start[i], &C_stop[i], &C[i], &trash1, &trash2)!=7) printf("Error1: %s\n", line);
    } else {i--;}
  }

  return 0;
}

/// A function to read the SDSSpix objects file and store the results in the appropriate arrays. 
int read_SDSSpix_file(FILE *objects, double *bin_number, double *solid_angle, double *ra, double *dec)
{
  long i;
  double pixel_area;
  char line[MAXCHARS];

  total_galaxies = 0.0;
  omega = 0.0;
  for (i=0; i<bins; i++) {
    if (fgets(line, sizeof(line), objects)==NULL) printf("Blank line. %s\n", line);
    if (line[0]!='#') {
      if (sscanf(line, "%lf %lf %lf %lf %lf", &bin_number[i], &solid_angle[i], &ra[i], &dec[i], &pixel_area)!=5) printf("Error1: %s\n", line);
      if (solid_angle[i] != 0.0) {
	solid_angle[i] *= pixel_area;
	omega += solid_angle[i];
	total_galaxies += bin_number[i];
      } else {i--;}
    } else {i--;}
  }

  return 0;
}

/// A function to read the healpix map and write it into the overdensity array.
int read_Healpix_file(double *overdensity, float *healpix, double *ra, double *dec, long nside)
{
  long i, j;
  double theta, phi;

  j = 0;
  omega = bins*129600.0/(M_PI*nside2npix(nside));
  for (i=0; i<nside2npix(nside); i++) {
    if (healpix[i] >= -1.0) {
      pix2ang_nest(nside, i, &theta, &phi);
      ra[j] = phi/d2r;
      dec[j] = 90.0-theta/d2r;
      if (ra[j] < 0.0) ra[j] += 360.0;
      overdensity[j] = healpix[i];
      j++;
    }
  }
  printf("#Area: %lf Pixels: %ld Total Galaxies = %ld\n", omega, bins, (long)total_galaxies);

  return 0;
}

/// A function to read the healpix map and write it into the overdensity array.
int read_cross_Healpix_file(double *overdensity, float *healpix_1, float *healpix_2, double *ra, double *dec, long nside)
{
  long i, j;
  double theta, phi;

  j = 0;
  omega = bins*129600.0/(M_PI*nside2npix(nside));
  for (i=0; i<nside2npix(nside); i++) {
    if ((healpix_1[i]>=-1.0)&&(healpix_2[i]>=-1.0)) {
      pix2ang_nest(nside, i, &theta, &phi);
      ra[j] = phi/d2r;
      dec[j] = 90.0-theta/d2r;
      if (ra[j] < 0.0) ra[j] += 360.0;
      overdensity[j] = healpix_1[i];
      j++;
    }
  }
  printf("#Area: %lf Pixels: %ld\n", omega, bins);

  return 0;
}

/** A function to calculate the signal and covariance matrices needed for the iterative process. 
  *	Returns the time taken to generate the signal matrix.  */
double calculate_matrices(double *bin_number, double *solid_angle, double *ra, double *dec, double *overdensity, double *cos_angle, double *noise, double *data_covariance, double *signal, int *C_start, int *C_stop)
{
  double signaltime;

  // Calculate overdensities and solid angle per pixel. 
  calculate_overdensity(overdensity, bin_number, solid_angle);

  // Calculate the covariance matrix based on the observed overdensities. 
  calculate_covariance(cos_angle, noise, data_covariance, ra, dec, solid_angle, overdensity);

  // Calculate the signal matrix using the angles between the pixels. 
  signaltime = calculate_signal(signal, cos_angle, C_start, C_stop);

  return signaltime;
}

/// A function to do KL compression on input data vectors and matrices.
int KL_compression(double *overdensity, double *signal, double *noise, double *data_covariance, double *C, FILE *output_KL)
{
  int i, j, l, INFO, LDA, LWORK, N;
  long k=0;
  double *A, *W, *WORK, KLtime = 0.0, *test1, *test2, *test3, *test4, *test5, *test6, *test7, *B, *Bt, *sum, *B_prime, *B_primet, total_snr = 0.0, signal_sum = 0.0;
  char JOBZ, UPLO;
  time_t t0, t1;

  time(&t0);

  N = LDA = bins;
  LWORK = 3*N - 1;
  JOBZ = 'V';  // 'N' for eigenvalues only, 'V' for eigenvalues and eigenvectors.
  UPLO = 'U';  // 'U' for upper triangle of B, 'L' for lower triangle of B.

  A = (double *)malloc(LDA*N*sizeof(double));
  B = (double *)malloc(LDA*N*sizeof(double));
  Bt = (double *)malloc(LDA*N*sizeof(double));
  W = (double *)malloc(N*sizeof(double));
  WORK = (double *)malloc(LWORK*sizeof(double)); 
  test1 = (double *)malloc(LDA*sizeof(double));
  test2 = (double *)malloc(LDA*N*sizeof(double));
  test3 = (double *)malloc(N*sizeof(double));
  test4 = (double *)malloc(LDA*N*sizeof(double));
  test5 = (double *)malloc(LDA*N*sizeof(double));
  test6 = (double *)malloc(LDA*N*sizeof(double));
  sum = (double *)malloc(bins*bins*sizeof(double));
 
  for (i=0; i<bins; i++) {
    for (j=0; j<bins; j++) {
      sum [i+j*bins]= 0.0;
      for (l=0; l<bands; l++) {
	sum[i+j*bins] += signal[i+j*bins+l*bins*bins]*C[l];
      }
    }
  }

  // A = inverse noise matrix times the signal matrix.  N^-1*S 
  invert_matrix(noise, bins); 
  multiply_matrices(noise, sum, B, bins, bins, bins);
  invert_matrix(noise, bins); 

  // Calculate eigenvalues and eigenvectors.
  printf("#Calculating  eigenvalues and eigenvectors.\n");
  dsyev(&JOBZ, &UPLO, &N, B, &LDA, W, WORK, &LWORK, &INFO);
  assert(INFO==0);

  // Reorder eigenvectors and eigenvalues from ascending to descending.
  for (i=0; i<bins; i++) {
    test3[i] = W[i];
    for (j=0; j<bins; j++) {
      test4[j+i*bins] = B[j+i*bins];
    }
  }
  for (i=0; i<bins; i++) {
    W[i] = test3[bins-1-i];
    for (j=0; j<bins; j++) {
      B[j+i*bins] = test4[j+(bins-1-i)*bins];
    }
    total_snr += W[i];
  }

  // Renormalizing so that Bt*N*B = 1.
  for (j=0; j<bins; j++) {
    for (i=0; i<bins; i++) {
      B[i+j*bins] = B[i+j*bins]/sqrt(noise[j+j*bins]);
    }
  }

  matrix_transpose(B, Bt, bins, bins);

  printf("#Printing out eigenvalues.\n");
  fprintf(output_KL, "#Printing out eigenvalues.\n");
  while (W[k] > 1.0) {
  //  printf("#Keeping %lf%% of the total signal to noise.\n", 100.0*SNR_RETAINED);
  //  while (signal_sum < total_snr*SNR_RETAINED) {
    printf("%ld %lf\n", k, W[k]);
    fprintf(output_KL, "%ld %lf\n", k, W[k]);
    signal_sum += W[k];
    k++;
  }

  B_prime = (double *)malloc(LDA*k*sizeof(double));
  B_primet = (double *)malloc(LDA*k*sizeof(double));
  test7 = (double *)malloc(k*k*bands*sizeof(double));

  // Create B_prime out of the leftmost columns of B, corresponding to eigenvalues greater than 1.
  for (i=0; i<k; i++) {
    for (j=0; j<bins; j++) {
      B_prime[j+i*bins] = B[j+i*bins];
    }
  }

  matrix_transpose(B_prime, B_primet, bins, k);

  // KL compress the data vector.
  printf("#Printing out new data vector y\n");
  fprintf(output_KL, "#Printing out new data vector y\n");
  multiply_matrices(B_primet, overdensity, test1, k, bins, 1);  
  for (i=0; i<bins; i++) {
    overdensity[i] = 0.0;
  }
  for (i=0; i<k; i++) {
    overdensity[i] = test1[i];
    printf("%d %lf\n", i, overdensity[i]);
    fprintf(output_KL, "%d %lf\n", i, overdensity[i]);
  }
  
  // KL compress the noise matrix.
  multiply_matrices(B_primet, noise, test5, k, bins, bins);  
  multiply_matrices(test5, B_prime, test6, k, bins, k);
  for (i=0; i<bins*bins; i++) {
    noise[i] = 0.0;
  }
  for (i=0; i<k*k; i++) {
    noise[i] = test6[i];
  }

  // KL compress the P matrices, denoted here by signal.
  for (l=0; l<bands; l++) {
    multiply_matrices(B_primet, &signal[l*bins*bins], test5, k, bins, bins);  
    multiply_matrices(test5, B_prime, &test7[l*k*k], k, bins, k);
  }

  for (i=0; i<bins*bins*bands; i++) {
    signal[i] = 0;
  }
  for (i=0; i<k*k*bands; i++) {
    signal[i] = test7[i];
  }

  // Reset bins to the number of modes with SNR > 1.0
  bins = k;
  printf("#New bins = %ld\n", k);

  // Recalculate the data_covariance matrix based on new data vector.
  for (j=0; j<bins; j++) {
    for (i=0; i<bins; i++) {
      data_covariance[i+j*bins] = overdensity[i]*overdensity[j]/((1.0-contamination)*(1.0-contamination));
    }
  } 
  printf("#Covariance matrix calculated from data.\n");

  time(&t1);
  KLtime = difftime(t1, t0);
  printf("#Done calculating KL-compression. Time = %lf seconds.  LWORK used: %d. Optimal LWORK: %lf.\n", KLtime, LWORK, WORK[0]);

  return 0;
}

/** A function to estimate the angular power spectrum given the signal and covariance matrices.
  * Returns the time taken by the inversion and matrix multiplication steps. */
double estimate_C(double *signal, double *model_covariance, double *data_covariance, double *noise, double *difference, double *average, double *A, double *B, double *F, double *C, int *C_start, int *C_stop, double *C_change, int iteration, FILE *output_C, FILE *output_Fisher, FILE *output_Window)
{
  double inversetime = 0.0, multiplicationtime = 0.0;

  printf("#Beginning recalculation of C[l] iteration number %d\n", iteration);
  
  // Calculate the expected covariance matrix based on C_l. 
  calculate_difference(signal, model_covariance, data_covariance, noise, difference, C);
  
  // Invert the model covariance matrix, result is saved in model_covariance. 
  inversetime += invert_matrix(model_covariance, bins);
  
  // Find the products of the model covariance matrix and derivative of the covariance matrix with respect to the bandpowers. 
  multiplicationtime += calculate_products(model_covariance, signal, A, B, average, difference);
  
  // Calculate the Fisher matrix and weighted average. 
  calculate_Fisher(F, A);
  
  // Print out the Fisher matrix. 
  print_Fisher(F, iteration, output_Fisher);
  
  // Invert the Fisher matrix. 
  invert_matrix(F, bands);
  
  // Recalculate the C_l. 
  calculate_KL_C(C, C_change, F, average, output_Window);
  
  // Print out the final values for this iteration. 
  print_values(C, C_start, C_stop, C_change, F, iteration, output_C);

  return (multiplicationtime+inversetime);
}

/// A function to calculate the overdensity from the pixelized galaxy counts and areas. 
int calculate_overdensity(double *overdensity, double *bin_number, double *solid_angle)
{
  long i;

  for (i=0; i<bins; i++) {
    overdensity[i] = bin_number[i]*omega/(total_galaxies*solid_angle[i])-1.0;
  }
  printf("#Overdensities calculated in %ld bins for %lf objects.\n", bins, total_galaxies);

  return 0;
}

/// A function to calulate the data and noise covariance matrices, as well as the cos_angle matrix. 
int calculate_covariance(double *cos_angle, double *noise, double *data_covariance, double *ra, double *dec, double *solid_angle, double *overdensity)
{
  long i, j;

  for (j=0; j<bins; j++) {
    for (i=0; i<bins; i++) {
      cos_angle[i+j*bins] = sin(dec[i]*d2r)*sin(dec[j]*d2r)+cos(dec[i]*d2r)*cos(dec[j]*d2r)*cos((ra[i]-ra[j])*d2r);
      noise[i+j*bins] = (i==j) ? (1.0*omega)/(1.0*total_galaxies)+LARGE_NUMBER : LARGE_NUMBER;
      data_covariance[i+j*bins] = overdensity[i]*overdensity[j]/((1.0-contamination)*(1.0-contamination));
    }
  } 
  printf("#Covariance matrix calculated from data.\n");

  return 0;
}

/// A function to calulate the data and noise covariance matrices, as well as the cos_angle matrix. 
int calculate_Healpix_covariance(double *cos_angle, double *noise, double *data_covariance, double *ra, double *dec, double *overdensity)
{
  long i, j;

  for (j=0; j<bins; j++) {
    for (i=0; i<bins; i++) {
      cos_angle[i+j*bins] = sin(dec[i]*d2r)*sin(dec[j]*d2r)+cos(dec[i]*d2r)*cos(dec[j]*d2r)*cos((ra[i]-ra[j])*d2r);
      noise[i+j*bins] = (i==j) ? (1.0*omega)/(1.0*total_galaxies)+LARGE_NUMBER : LARGE_NUMBER;
      data_covariance[i+j*bins] = overdensity[i]*overdensity[j]/((1.0-contamination)*(1.0-contamination));
    }
  } 
  printf("#Covariance matrix calculated from data.\n");

  return 0;
}

/// A function to calulate the data and noise covariance matrices, as well as the cos_angle matrix. 
int calculate_cross_covariance(double *cos_angle, double *noise, double *data_covariance, double *ra, double *dec, double *overdensity1, double *overdensity2)
{
  long i, j;

  for (j=0; j<bins; j++) {
    for (i=0; i<bins; i++) {
      cos_angle[i+j*bins] = sin(dec[i]*d2r)*sin(dec[j]*d2r)+cos(dec[i]*d2r)*cos(dec[j]*d2r)*cos((ra[i]-ra[j])*d2r);
      noise[i+j*bins] = 0.0;
      data_covariance[i+j*bins] = overdensity1[i]*overdensity2[j]/((1.0-contamination)*(1.0-contamination));
    }
  } 
  printf("#Cross covariance matrix calculated from data.\n");

  return 0;
}

/// A function to calulate the cos_angle matrix. 
int calculate_cos_angle(double *cos_angle, double *ra, double *dec)
{
  long i, j;

  for (j=0; j<bins; j++) {
    for (i=0; i<bins; i++) {
      cos_angle[i+j*bins] = sin(dec[i]*d2r)*sin(dec[j]*d2r)+cos(dec[i]*d2r)*cos(dec[j]*d2r)*cos((ra[i]-ra[j])*d2r);
    }
  } 
  printf("#Cos_angle matrix calculated from data.\n");

  return 0;
}

///  A parallelized function to calculate the signal matrix using the cos_angle matrix and bandpower limits. 
double calculate_signal(double *signal, double *cos_angle, int *C_start, int *C_stop)
{
  long i, j, k, l, z;
  double signaltime;
  time_t t0, t1;

  printf("#Calculating signal matrix.\n");
  time(&t0);
#pragma omp parallel for shared(bands, bins, signal, C_start, C_stop, cos_angle) private(l, j, i, k)
  for (z=0; z<num_proc; z++) {
    for (i=z; i<(bands*bins*bins); i+=num_proc) {
      j = i%(bins*bins);
      l = i/(bins*bins);
      signal[i] = 0.0;
      for (k=C_start[l]; k<=C_stop[l]; k++) {
	signal[i] += ((2.0*k+1.0)/(2.0*k*(k+1)))*Legendre(cos_angle[j], k);
      }
    }
  }
  time(&t1);
  signaltime = difftime(t1, t0);
  printf("#Done calculating signal matrix. Time = %lf seconds.\n", signaltime);

  return signaltime;
}

int print_signal(FILE *output_signal, double *signal)
{
  long i;

  for (i=0; i<(bands*bins*bins); i++) {
    fprintf(output_signal, "%lf\n", signal[i]);
  }
  printf("#Done printing signal matrix.\n");  

  return 0;
}

int read_signal(FILE *input_signal, double *signal)
{
  long i;
  char line[MAXCHARS];
  time_t t0, t1;

  printf("#Reading signal matrix.\n");
  time(&t0);
  for (i=0; i<(bands*bins*bins); i++){
    if (fgets(line, sizeof(line), input_signal)==NULL) printf("Blank line. %s\n", line);
    if (sscanf(line, "%lf", &signal[i])!=1) printf("Error1: %s\n", line);
  }
  time(&t1);
  printf("#Done calculating signal matrix. Time = %lf seconds.\n", difftime(t1, t0));

  return 0;
}

/// A function to calculate the difference of the data and model covariance matrices. 
int calculate_difference(double *signal, double *model_covariance, double *data_covariance, double *noise, double *difference, double *C)
{
  long i, j, l;
  double sum;

  for (j=0; j<bins; j++) {
    for (i=0; i<bins; i++) {
      sum = 0.0;
      for (l=0; l<bands; l++) {
	sum += signal[i+j*bins+l*bins*bins]*C[l];
      }
      model_covariance[i+j*bins] = sum+noise[i+j*bins];
      difference[i+j*bins] = data_covariance[i+j*bins]-noise[i+j*bins];
    }
  }
  
  return 0;
}

/**  A function to calculate the product of the derivative of the signal matrix with the model covariance matrix.
  *	Also calculates the product of the above result with another factor of the model covariance and the traces of the 
  *	multiplications of the first results. Returns the time it takes to do these multiplications. */
double calculate_products(double *model_covariance, double *signal, double *A, double *B, double *average, double *difference)
{
  int l, LDA, LDB, LDC, M, N, K;
  double ALPHA, BETA, multiplicationtime;
  char TRANSA, TRANSB;
  time_t t0, t1;

  M = N = K = LDA = LDB = LDC = bins;
  ALPHA = 1.0;
  BETA = 0.0;
  TRANSA = TRANSB = 'N';

  time(&t0);
  for (l=0; l<bands; l++) {
    dgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, model_covariance, &LDA, &signal[l*bins*bins], &LDB, &BETA, &A[l*bins*bins], &LDC);
    dgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, &A[l*bins*bins], &LDA, model_covariance, &LDB, &BETA, &B[l*bins*bins], &LDC);
    average[l] = trace_multiply(difference, &B[l*bins*bins], bins, bins);
  }
  time(&t1);
  multiplicationtime = difftime(t1, t0);
  printf("#Done calculating matrix multiplications. Time = %lf seconds.\n", multiplicationtime);    
  
  return multiplicationtime;
}

/// A function to calculate the Fisher matrix. 
int calculate_Fisher(double *F, double *A)
{
  int l, l_prime;

  for (l=0; l<bands; l++) {
    for (l_prime=0; l_prime<=l; l_prime++) {
      F[l+l_prime*bands] = 0.5*trace_multiply(&A[l*bins*bins], &A[l_prime*bins*bins], bins, bins);
      F[l_prime+l*bands] = F[l+l_prime*bands];
    }
  }

  return 0;
}

/// A function to print the Fisher matrix to stdout. 
int print_Fisher(double *F, int iteration, FILE *output_Fisher)
{
  int l, l_prime;

  printf("#Printing out Fisher Information matrix for iteration %d.\n", iteration);
  for (l=0; l<bands; l++) {
    for (l_prime=0; l_prime<bands; l_prime++) {
      printf("%e ", F[l+l_prime*bands]);
      fprintf(output_Fisher, "%e ", F[l+l_prime*bands]);
    }
    printf("\n");
    fprintf(output_Fisher, "\n");
  }
  fflush(NULL);

  return 0;
}

/** A function to calculate the C_l from the Fisher matrix. Returns non-zero 
  *	if the change in C_l are greater than the errors. */
int calculate_C(double *C, double *C_change, double *F, double *average)
{
  int l, l_prime, hold;

  hold = 0;
  for (l=0; l<bands; l++) {
    C_change[l] = 0.0;
    for (l_prime = 0; l_prime<bands; l_prime++) {
      C_change[l] += 0.5*F[l+l_prime*bands]*average[l_prime];
    }
    if (fabs(C_change[l]) > fabs(sqrt(F[l+l*bands]))) hold += 1;
    C[l] = C_change[l];
  }

  return hold;
}

/** A function to calculate the C_l from the Fisher matrix. Returns non-zero 
  * if the change in C_l are greater than the errors. */
int calculate_KL_C(double *C, double *C_change, double *F, double *average, FILE *output_Window)
{
  int l, l_prime, hold;
  double sum, *invsqrtF, *A, *B, *D, *W;

  invsqrtF = (double *)malloc(bands*bands*sizeof(double));
  A = (double *)malloc(bands*bands*sizeof(double));
  B = (double *)malloc(bands*bands*sizeof(double));
  D = (double *)malloc(bands*bands*sizeof(double));
  W = (double *)malloc(bands*bands*sizeof(double));

  matrix_square_root(F, invsqrtF, bands);
  invert_matrix(F, bands);
  multiply_matrices(invsqrtF, F, A, bands, bands, bands); // A = F^1/2

  for (l_prime=0; l_prime<bands; l_prime++) {
    sum = 0.0;
    for (l=0; l<bands; l++) {
      sum += A[l_prime+l*bands];
    }
    for (l=0; l<bands; l++) {
      W[l_prime+l*bands] = A[l_prime+l*bands]/sum;
    }  
  }

  invert_matrix(A, bands);
  multiply_matrices(W, A, D, bands, bands, bands); // D = WF^-1/2
  invert_matrix(A, bands);
  multiply_matrices(D, invsqrtF, B, bands, bands, bands);
  multiply_matrices(B, F, W, bands, bands, bands);
  invert_matrix(F, bands);

  printf("#Printing out Window matrix.\n");  
  for (l_prime=0; l_prime<bands; l_prime++) {
    for (l=0; l<bands; l++) {
      printf("%lf ", W[l*bands+l_prime]);
      fprintf(output_Window, "%lf ", W[l*bands+l_prime]);
    }
    printf("\n");
    fprintf(output_Window, "\n");
  }

  hold = 0;
  for (l=0; l<bands; l++) {
    C_change[l] = 0.0;
    for (l_prime = 0; l_prime<bands; l_prime++) {
      C_change[l] += 0.5*B[l+l_prime*bands]*average[l_prime];
    }
    if (fabs(C_change[l]) > fabs(sqrt(F[l+l*bands]))) hold += 1;
    C[l] = C_change[l];
  }

  return hold;
}

/// A function to print the C_l values and errors to stdout. 
int print_values(double *C, int *C_start, int *C_stop, double *C_change, double *F, int iteration, FILE *output_C)
{
  int l;

  printf("#All C[l] have been recalculated for iteration = %d.\n", iteration);
  for (l=0; l<bands; l++) {
    printf("%d %d %d %d %lf %lf %lf\n", l, (C_stop[l]+C_start[l])/2, C_start[l], C_stop[l], C[l], sqrt(F[l+l*bands]), C_change[l]);
    fprintf(output_C, "%d %d %d %d %lf %lf %lf\n", l, (C_stop[l]+C_start[l])/2, C_start[l], C_stop[l], C[l], sqrt(F[l+l*bands]), C_change[l]);
  }

  return 0;
}
