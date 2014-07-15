#include "angular_power_spectrum.h"
//#include "new_angular_power_spectrum.h"

/// Open the Healpix file and count the pixels and galaxies.
long count_Healpix_pixels(char *input_healpix, float *healpix)
{
  int i;
  long NSIDE;
  char *trash[kMaxChars], name[kMaxChars]; // ordering[10], coords[1];

  // Strip the directory structure from the filename if present
  if (strrchr(input_healpix, '/')!=NULL) {
    strcpy(name, strrchr(input_healpix, '/')+1);
  } else {
    strcpy(name, input_healpix);
  }

  NSIDE = strtol(strpbrk(name, "0123456789"), trash, 10);
  g_total_galaxies = strtod(strpbrk(*trash, "0123456789"), (char**)NULL);
  printf("#Done reading %ld lines.\n", nside2npix(NSIDE));

  g_bins = 0;
  for (i=0; i<nside2npix(NSIDE); i++) {
    if (healpix[i] >= -1.0) {
      g_bins++;
    }
  }
  printf("#Done counting %ld pixels.\n", g_bins);

  return NSIDE;
}

/// A function to read the bandpowers file, and store the results in the appropriate arrays. 
int read_bandpower_file(FILE *bandpowers, double *C, int *C_start, int *C_stop)
{
  int i, trash, garbage;
  double trash1, trash2;
  char line[kMaxChars];

  for (i=0; i<g_bands; i++) {
    if (fgets(line, sizeof(line), bandpowers)==NULL) printf("Blank line. %s\n", line);
    if (line[0]!='#') {
      if (sscanf(line, "%d %d %d %d %lf %lf %lf", &trash, &garbage, &C_start[i], &C_stop[i], &C[i], &trash1, &trash2)!=7) printf("Error1: %s\n", line);
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
  g_omega = g_bins*129600.0/(M_PI*nside2npix(nside));
  for (i=0; i<nside2npix(nside); i++) {
    if (healpix[i] >= -1.0) {
      pix2ang_nest(nside, i, &theta, &phi);
      ra[j] = phi/kDegreeToRadian;
      dec[j] = 90.0-theta/kDegreeToRadian;
      if (ra[j] < 0.0) ra[j] += 360.0;
      overdensity[j] = healpix[i];
      j++;
    }
  }
  printf("#Area: %lf Pixels: %ld Total Galaxies = %ld\n", g_omega, g_bins, (long)g_total_galaxies);

  return 0;
}


/// A function to do KL compression on input data vectors and matrices.
int KL_compression(double *overdensity, double *signal, double *noise, double *data_covariance, double *C, FILE *output_KL, char* test_root)
{
  int i, j, l, INFO, LDA, LWORK, N;
  long k=0;
  double *A, *W, *WORK, KLtime = 0.0, *test1, *test2, *test3, *test4, *test5, *test6, *test7, *B, *Bt, *sum, *B_prime, *B_primet, total_snr = 0.0, signal_sum = 0.0;
  char JOBZ, UPLO;
  Timer compression_time;
  tic(&compression_time);

  N = LDA = g_bins;
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
  sum = (double *)malloc(g_bins*g_bins*sizeof(double));
 
  for (i=0; i<g_bins; i++) {
    for (j=0; j<g_bins; j++) {
      sum [i+j*g_bins]= 0.0;
      for (l=0; l<g_bands; l++) {
        sum[i+j*g_bins] += signal[i+j*g_bins+l*g_bins*g_bins]*C[l];
      }
    }
  }

  // A = inverse noise matrix times the signal matrix.  N^-1*S 
  // printf("noise: %g %g\n", noise[0]-KLargeNumber, noise[1]-KLargeNumber);
  // for (i = 0; i < g_bins; ++i){
  //   for (j = 0; j < g_bins; ++j){
  //     printf("%g ", noise[j+i*g_bins]);
  //   }
  //   printf("\n");
  // }
  invert_matrix(noise, g_bins);
  printf("inverse noise: %f %f\n", noise[0], noise[1]);
  multiply_matrices(noise, sum, B, g_bins, g_bins, g_bins);
  invert_matrix(noise, g_bins);

  // //Test sym
  // printf("#-#########\n#-#Testing symmetry\n#-#########\n");
  // double threshold = .01;
  // double diff;
  // int count = 0;
  // for (i=0; i<g_bins; i++) {
  //   for (j=i+1; j<g_bins; j++) {
  //     diff = B[i+j*g_bins] - B[j+i*g_bins];
  //     diff = (diff < 0)? -diff : diff;
  //     if (diff>threshold) {
  //       printf("#-# %10f =/= %9f    diff: %8f    diff/avg: %9f\n", B[i+j*g_bins], B[j+i*g_bins], diff, 2.0*diff/(B[i+j*g_bins]+B[j+i*g_bins]) );
  //       ++count;
  //     }
  //     if (count > 100) break;
  //   }
  //   if (count > 100) break;
  // }
  // printf("#-# With threshold %g\n", threshold);
  // printf("#-# %d out of %.0f don't match\n", count, 0.5*(g_bins*g_bins-g_bins));
  // printf("#-# %.2f percent\n", 100.0*count/(0.5*(g_bins*g_bins-g_bins)));
  // printf("#-#########\n#-#Finished Testing\n#-#########\n");



  // Calculate eigenvalues and eigenvectors.
  printf("#Calculating  eigenvalues and eigenvectors.\n");
  dsyev(&JOBZ, &UPLO, &N, B, &LDA, W, WORK, &LWORK, &INFO);
  assert(INFO==0);

  // Reorder eigenvectors and eigenvalues from ascending to descending.
  for (i=0; i<g_bins; i++) {
    test3[i] = W[i];
    for (j=0; j<g_bins; j++) {
      test4[j+i*g_bins] = B[j+i*g_bins];
    }
  }
  for (i=0; i<g_bins; i++) {
    W[i] = test3[g_bins-1-i];
    for (j=0; j<g_bins; j++) {
      B[j+i*g_bins] = test4[j+(g_bins-1-i)*g_bins];
    }
    total_snr += W[i];
  }

  // Renormalizing so that Bt*N*B = 1.
  for (j=0; j<g_bins; j++) {
    for (i=0; i<g_bins; i++) {
      B[i+j*g_bins] = B[i+j*g_bins]/sqrt(noise[j+j*g_bins]);
    }
  }

# ifdef APS_OUTPUT_TEST
  save_raw_double_array(test_root, "eigenvectors", B, g_bins * g_bins);
# endif

  matrix_transpose(B, Bt, g_bins, g_bins);

  printf("#Printing out eigenvalues.\n");
  fprintf(output_KL, "#Printing out eigenvalues.\n");
  while (W[k] > 1.0) {
  //  printf("#Keeping %lf%% of the total signal to noise.\n", 100.0*kSnrRetained);
  //  while (signal_sum < total_snr*kSnrRetained) {
#   ifndef APS_SUPPRESS_MATRIX_STDOUT
    printf("%ld %lf\n", k, W[k]);
#   endif
    fprintf(output_KL, "%ld %lf\n", k, W[k]);
    signal_sum += W[k];
    k++;
  }

  B_prime = (double *)malloc(LDA*k*sizeof(double));
  B_primet = (double *)malloc(LDA*k*sizeof(double));
  test7 = (double *)malloc(k*k*g_bands*sizeof(double));

  // Create B_prime out of the leftmost columns of B, corresponding to eigenvalues greater than 1.
  // Row major is: row*NUMCOLS + col | i is column, so this is in Column major
  for (i=0; i<k; i++) {
    for (j=0; j<g_bins; j++) {
      B_prime[j+i*g_bins] = B[j+i*g_bins];
    }
  }

  matrix_transpose(B_prime, B_primet, g_bins, k);

  // KL compress the data vector.
  printf("#Printing out new data vector y\n");
  fprintf(output_KL, "#Printing out new data vector y\n");
  multiply_matrices(B_primet, overdensity, test1, k, g_bins, 1);  
  for (i=0; i<g_bins; i++) {
    overdensity[i] = 0.0;
  }
  for (i=0; i<k; i++) {
    overdensity[i] = test1[i];
#   ifndef APS_SUPPRESS_MATRIX_STDOUT
    printf("%d %lf\n", i, overdensity[i]);
#   endif
    fprintf(output_KL, "%d %lf\n", i, overdensity[i]);
  }
  
  // KL compress the noise matrix.
  multiply_matrices(B_primet, noise, test5, k, g_bins, g_bins);  
  multiply_matrices(test5, B_prime, test6, k, g_bins, k);
  for (i=0; i<g_bins*g_bins; i++) {
    noise[i] = 0.0;
  }
  for (i=0; i<k*k; i++) {
    noise[i] = test6[i];
  }

  // KL compress the P matrices, denoted here by signal.
  for (l=0; l<g_bands; l++) {
    multiply_matrices(B_primet, &signal[l*g_bins*g_bins], test5, k, g_bins, g_bins);  
    multiply_matrices(test5, B_prime, &test7[l*k*k], k, g_bins, k);
  }

  for (i=0; i<g_bins*g_bins*g_bands; i++) {
    signal[i] = 0;
  }
  for (i=0; i<k*k*g_bands; i++) {
    signal[i] = test7[i];
  }

  // Reset g_bins to the number of modes with SNR > 1.0
  g_bins = k;
  printf("#New g_bins = %ld\n", k);

  // Recalculate the data_covariance matrix based on new data vector.
  for (j=0; j<g_bins; j++) {
    for (i=0; i<g_bins; i++) {
      data_covariance[i+j*g_bins] = overdensity[i]*overdensity[j]/((1.0-kContamination)*(1.0-kContamination));
    }
  } 
  printf("#Covariance matrix calculated from data.\n");


  toc(&compression_time);
  printf("#Done calculating KL-compression. Time = %g seconds.  LWORK used: %d. Optimal LWORK: %g.\n", compression_time.elapsed, LWORK, WORK[0]);

  return 0;
}

/** A function to estimate the angular power spectrum given the signal and covariance matrices.
  * Returns the time taken by the inversion and matrix multiplication steps. */
double estimate_C(double *signal, double *model_covariance, double *data_covariance, double *noise, double *difference, double *average, double *A, double *B, double *F, double *C, int *C_start, int *C_stop, double *C_change, int iteration, FILE *output_C, FILE *output_Fisher, FILE *output_Window, char* test_root)
{
  double inversetime = 0.0, multiplicationtime = 0.0;

  printf("#Beginning recalculation of C[l] iteration number %d\n", iteration);
  
  // Calculate the expected covariance matrix based on C_l. 
  calculate_difference(signal, model_covariance, data_covariance, noise, difference, C);
  
  // Invert the model covariance matrix, result is saved in model_covariance. 
  inversetime += invert_matrix(model_covariance, g_bins);
  
  // Find the products of the model covariance matrix and derivative of the covariance matrix with respect to the bandpowers. 
  multiplicationtime += calculate_products(model_covariance, signal, A, B, average, difference);
  
  // Calculate the Fisher matrix and weighted average. 
  calculate_Fisher(F, A);
# ifdef APS_OUTPUT_TEST
  char filename[kMaxChars];
  sprintf(filename, "iter_%d_fisher", iteration);
  save_raw_double_array(test_root, filename, F, g_bands*g_bands);
# endif
  
  // Print out the Fisher matrix. 
  print_Fisher(F, iteration, output_Fisher);
  
  // Invert the Fisher matrix. 
  invert_matrix(F, g_bands);
  
  // Recalculate the C_l. 
  calculate_KL_C(C, C_change, F, average, output_Window);
  
  // Print out the final values for this iteration. 
  print_values(C, C_start, C_stop, C_change, F, iteration, output_C);

  return (multiplicationtime+inversetime);
}


/// A function to calulate the data and noise covariance matrices, as well as the cos_angle matrix. 
int calculate_Healpix_covariance(double *cos_angle, double *noise, double *data_covariance, double *ra, double *dec, double *overdensity)
{
  long i, j;

  for (j=0; j<g_bins; j++) {
    for (i=0; i<g_bins; i++) {
      cos_angle[i+j*g_bins] = sin(dec[i]*kDegreeToRadian)*sin(dec[j]*kDegreeToRadian)+cos(dec[i]*kDegreeToRadian)*cos(dec[j]*kDegreeToRadian)*cos((ra[i]-ra[j])*kDegreeToRadian);
      noise[i+j*g_bins] = (i==j) ? (1.0*g_omega)/(1.0*g_total_galaxies)+KLargeNumber : KLargeNumber;
      data_covariance[i+j*g_bins] = overdensity[i]*overdensity[j]/((1.0-kContamination)*(1.0-kContamination));
    }
  } 
  printf("#Covariance matrix calculated from data.\n");

  return 0;
}


///  A parallelized function to calculate the signal matrix using the cos_angle matrix and bandpower limits. 
double calculate_signal(double *signal, double *cos_angle, int *C_start, int *C_stop)
{
  long i, j, k, l, z;
  Timer time_calc_signal;
  tic(&time_calc_signal);

  printf("#Calculating signal matrix.\n");

#pragma omp parallel for shared(g_bands, g_bins, signal, C_start, C_stop, cos_angle) private(l, j, i, k)
  for (z=0; z<kNumThreads; z++) {
    for (i=z; i<(g_bands*g_bins*g_bins); i+=kNumThreads) {
      j = i%(g_bins*g_bins);
      l = i/(g_bins*g_bins);
      signal[i] = 0.0;
      for (k=C_start[l]; k<=C_stop[l]; k++) {
        signal[i] += ((2.0*k+1.0)/(2.0*k*(k+1)))*Legendre(cos_angle[j], k);
      }
    }
  }
  toc(&time_calc_signal);
  printf("#Done calculating signal matrix. Time = %g seconds.\n", time_calc_signal.elapsed);

  return time_calc_signal.elapsed;
}

int print_signal(FILE *output_signal, double *signal)
{
  long i;

  for (i=0; i<(g_bands*g_bins*g_bins); i++) {
    fprintf(output_signal, "%lf\n", signal[i]);
  }
  printf("#Done printing signal matrix.\n");  

  return 0;
}

/// A function to calculate the difference of the data and model covariance matrices. 
int calculate_difference(double *signal, double *model_covariance, double *data_covariance, double *noise, double *difference, double *C)
{
  long i, j, l;
  double sum;

  for (j=0; j<g_bins; j++) {
    for (i=0; i<g_bins; i++) {
      sum = 0.0;
      for (l=0; l<g_bands; l++) {
        sum += signal[i+j*g_bins+l*g_bins*g_bins]*C[l];
      }
      model_covariance[i+j*g_bins] = sum+noise[i+j*g_bins];
      difference[i+j*g_bins] = data_covariance[i+j*g_bins]-noise[i+j*g_bins];
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
  Timer time_calc_products;
  tic(&time_calc_products);

  M = N = K = LDA = LDB = LDC = g_bins;
  ALPHA = 1.0;
  BETA = 0.0;
  TRANSA = TRANSB = 'N';

  
  for (l=0; l<g_bands; l++) {
    //A[l*g_bins*g_bins] := model_covariance * signal[l*g_bins*g_bins]
    dgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, model_covariance, &LDA, &signal[l*g_bins*g_bins], &LDB, &BETA, &A[l*g_bins*g_bins], &LDC);
    //B[l*g_bins*g_bins] := A[l*g_bins*g_bins] * model_covariance
    //B[l*g_bins*g_bins] := model_covariance * signal[l*g_bins*g_bins] * model_covariance
    dgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, &A[l*g_bins*g_bins], &LDA, model_covariance, &LDB, &BETA, &B[l*g_bins*g_bins], &LDC);
    average[l] = trace_multiply(difference, &B[l*g_bins*g_bins], g_bins, g_bins);
  }
  
  toc(&time_calc_products);
  printf("#Done calculating matrix multiplications. Time = %g seconds.\n", time_calc_products.elapsed);
  
  return time_calc_products.elapsed;
}

/// A function to calculate the Fisher matrix. 
int calculate_Fisher(double *F, double *A)
{
  int l, l_prime;

  for (l=0; l<g_bands; l++) {
    for (l_prime=0; l_prime<=l; l_prime++) {
      // F := trace(A[l] * A[l_prime])/2
      F[l+l_prime*g_bands] = 0.5*trace_multiply(&A[l*g_bins*g_bins], &A[l_prime*g_bins*g_bins], g_bins, g_bins);
      F[l_prime+l*g_bands] = F[l+l_prime*g_bands];
    }
  }

  return 0;
}

/// A function to print the Fisher matrix to stdout. 
int print_Fisher(double *F, int iteration, FILE *output_Fisher)
{
  int l, l_prime;

  printf("#Printing out Fisher Information matrix for iteration %d.\n", iteration);
  for (l=0; l<g_bands; l++) {
    for (l_prime=0; l_prime<g_bands; l_prime++) {
#     ifndef APS_SUPPRESS_MATRIX_STDOUT
      printf("%e ", F[l+l_prime*g_bands]);
#     endif
      fprintf(output_Fisher, "%e ", F[l+l_prime*g_bands]);
    }
    printf("\n");
    fprintf(output_Fisher, "\n");
  }
  fflush(NULL);

  return 0;
}


/** A function to calculate the C_l from the Fisher matrix. Returns non-zero 
  * if the change in C_l are greater than the errors. */
int calculate_KL_C(double *C, double *C_change, double *F, double *average, FILE *output_Window)
{
  int l, l_prime, hold;
  double sum, *invsqrtF, *A, *B, *D, *W;

  invsqrtF = (double *)malloc(g_bands*g_bands*sizeof(double));
  A = (double *)malloc(g_bands*g_bands*sizeof(double));
  B = (double *)malloc(g_bands*g_bands*sizeof(double));
  D = (double *)malloc(g_bands*g_bands*sizeof(double));
  W = (double *)malloc(g_bands*g_bands*sizeof(double));

  matrix_square_root(F, invsqrtF, g_bands);
  invert_matrix(F, g_bands);
  multiply_matrices(invsqrtF, F, A, g_bands, g_bands, g_bands); // A = F^1/2

  for (l_prime=0; l_prime<g_bands; l_prime++) {
    sum = 0.0;
    for (l=0; l<g_bands; l++) {
      sum += A[l_prime+l*g_bands];
    }
    for (l=0; l<g_bands; l++) {
      W[l_prime+l*g_bands] = A[l_prime+l*g_bands]/sum;
    }  
  }

  invert_matrix(A, g_bands);
  multiply_matrices(W, A, D, g_bands, g_bands, g_bands); // D = WF^-1/2
  invert_matrix(A, g_bands);
  multiply_matrices(D, invsqrtF, B, g_bands, g_bands, g_bands);
  multiply_matrices(B, F, W, g_bands, g_bands, g_bands);
  invert_matrix(F, g_bands);

  printf("#Printing out Window matrix.\n");  
  for (l_prime=0; l_prime<g_bands; l_prime++) {
    for (l=0; l<g_bands; l++) {
#     ifndef APS_SUPPRESS_MATRIX_STDOUT
      printf("%lf ", W[l*g_bands+l_prime]);
#     endif
      fprintf(output_Window, "%lf ", W[l*g_bands+l_prime]);
    }
    printf("\n");
    fprintf(output_Window, "\n");
  }

  hold = 0;
  for (l=0; l<g_bands; l++) {
    C_change[l] = 0.0;
    for (l_prime = 0; l_prime<g_bands; l_prime++) {
      C_change[l] += 0.5*B[l+l_prime*g_bands]*average[l_prime];
    }
    if (fabs(C_change[l]) > fabs(sqrt(F[l+l*g_bands]))) hold += 1;
    C[l] = C_change[l];
  }

  return hold;
}

/// A function to print the C_l values and errors to stdout. 
int print_values(double *C, int *C_start, int *C_stop, double *C_change, double *F, int iteration, FILE *output_C)
{
  int l;

  printf("#All C[l] have been recalculated for iteration = %d.\n", iteration);
  for (l=0; l<g_bands; l++) {
#   ifndef APS_SUPPRESS_MATRIX_STDOUT
    printf("%d %d %d %d %lf %lf %lf\n", l, (C_stop[l]+C_start[l])/2, C_start[l], C_stop[l], C[l], sqrt(F[l+l*g_bands]), C_change[l]);
#   endif
    fprintf(output_C, "%d %d %d %d %lf %lf %lf\n", l, (C_stop[l]+C_start[l])/2, C_start[l], C_stop[l], C[l], sqrt(F[l+l*g_bands]), C_change[l]);
  }

  return 0;
}