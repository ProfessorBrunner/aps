// Functions that are not used in the distributed
// APS code.


//
//  From angular_power_spectrum.c
///////////////////////////////////////////////////////////

/// A function to count the pixels with non-zero area in the objects file.  
int count_usable_pixels(FILE *objects, int number_lines)
{
  int i;
  double lam, et, trash1, trash2, pixel_area;
  char line[kMaxChars];

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
long count_cross_Healpix_pixels(char *input_healpix, float *healpix_1, float *healpix_2)
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

/// A function to read the SDSSpix objects file and store the results in the appropriate arrays. 
int read_SDSSpix_file(FILE *objects, double *bin_number, double *solid_angle, double *ra, double *dec)
{
  long i;
  double pixel_area;
  char line[kMaxChars];

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
int read_cross_Healpix_file(double *overdensity, float *healpix_1, float *healpix_2, double *ra, double *dec, long nside)
{
  long i, j;
  double theta, phi;

  j = 0;
  omega = bins*129600.0/(M_PI*nside2npix(nside));
  for (i=0; i<nside2npix(nside); i++) {
    if ((healpix_1[i]>=-1.0)&&(healpix_2[i]>=-1.0)) {
      pix2ang_nest(nside, i, &theta, &phi);
      ra[j] = phi/kDegreeToRadian;
      dec[j] = 90.0-theta/kDegreeToRadian;
      if (ra[j] < 0.0) ra[j] += 360.0;
      overdensity[j] = healpix_1[i];
      j++;
    }
  }
  printf("#Area: %lf Pixels: %ld\n", omega, bins);

  return 0;
}


/** A function to calculate the signal and covariance matrices needed for the iterative process. 
  * Returns the time taken to generate the signal matrix.  */
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
      cos_angle[i+j*bins] = sin(dec[i]*kDegreeToRadian)*sin(dec[j]*kDegreeToRadian)+cos(dec[i]*kDegreeToRadian)*cos(dec[j]*kDegreeToRadian)*cos((ra[i]-ra[j])*kDegreeToRadian);
      noise[i+j*bins] = (i==j) ? (1.0*omega)/(1.0*total_galaxies)+KLargeNumber : KLargeNumber;
      data_covariance[i+j*bins] = overdensity[i]*overdensity[j]/((1.0-kContamination)*(1.0-kContamination));
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
      cos_angle[i+j*bins] = sin(dec[i]*kDegreeToRadian)*sin(dec[j]*kDegreeToRadian)+cos(dec[i]*kDegreeToRadian)*cos(dec[j]*kDegreeToRadian)*cos((ra[i]-ra[j])*kDegreeToRadian);
      noise[i+j*bins] = 0.0;
      data_covariance[i+j*bins] = overdensity1[i]*overdensity2[j]/((1.0-kContamination)*(1.0-kContamination));
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
      cos_angle[i+j*bins] = sin(dec[i]*kDegreeToRadian)*sin(dec[j]*kDegreeToRadian)+cos(dec[i]*kDegreeToRadian)*cos(dec[j]*kDegreeToRadian)*cos((ra[i]-ra[j])*kDegreeToRadian);
    }
  } 
  printf("#Cos_angle matrix calculated from data.\n");

  return 0;
}


int read_signal(FILE *input_signal, double *signal)
{
  long i;
  char line[kMaxChars];
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

/** A function to calculate the C_l from the Fisher matrix. Returns non-zero 
  * if the change in C_l are greater than the errors. */
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

//
//  From math.c
///////////////////////////////////////////////////////////

/// Calculate A * B = C where A and B are diagonal square N x N matrices.
int multiply_diagonal_matrices(double *A, double *B, double *C, int N)
{
  int i;

  for (i=0; i<N; i++) {
    C[i+i*N] = A[i+i*N]*B[i+i*N];
  }

  return 0;
}


/// Calculate the trace of N x N matrix A
double matrix_trace(double *A, int N)
{
  int i;
  double trace;

  trace = 0.0;
  for (i=0; i<N; i++) {
    trace += A[i+i*N];
  }

  return trace;
}

/// Multiplication of two square size by size matrices: AB = C
double multiply_square_matrices(double *A, double *B, double *C, long size)
{
  int LDA, LDB, LDC, M, N, K;
  double ALPHA, BETA, multiplicationtime;
  char TRANSA, TRANSB;
  time_t t0, t1;

  M = N = K = LDA = LDB = LDC = size;
  ALPHA = 1.0;
  BETA = 0.0;
  TRANSA = TRANSB = 'N';

  time(&t0);
  dgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
  time(&t1);
  multiplicationtime = difftime(t1, t0);
  printf("#Done calculating matrix multiplication. Time = %lf seconds.\n", multiplicationtime);

  return multiplicationtime;
}


double dot_product(double *A, double *B, long size)
{
  long i;
  double total=0.0;

  for (i=0; i<size; i++) {
    total += A[i]*B[i];
  }

  return total;
}
