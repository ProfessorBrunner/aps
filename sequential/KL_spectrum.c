#include "angular_power_spectrum.h"
//#include "new_angular_power_spectrum.h"

// Run with "./KL_spectrum <data.dat> <bandpowers.dat>

/*
  Runtimes:
  ~700 pixels: 15 min walltime / 1 proc / ?? mem  (interactive)
  ~2700 pixels: 15 min walltime / 24 proc / ~4G mem
*/

long g_bins;               // Total number of usable pixels.
long g_bands;              // Total number of bandpowers.
double g_omega;            // Total area of usable pixels.
double g_total_galaxies;   // Total number of galaxies in usable pixels.

/** Given the pixelized galaxy file and bandpower file, this program will calculate the angular power 
	* spectrum that best fits the input files. */
int main(int argc, char *argv[])
{
  int i, n, *C_start, *C_stop;
  long nside, NSIDE;
  float *healpix;
  double *overdensity, *data_covariance, *noise, *cos_angle, *difference, *A, *B, *C, *C_change; 
  double *ra, *dec, *F, *average, *model_covariance, *signal, signaltime, iterationtime = 0;
  char ordering[10], coords[1], filename[kMaxChars], name[kMaxChars];
  time_t t0, t1, t2, t3;
  FILE *objects, *bandpowers, *output_KL, *output_C, *output_Fisher, *output_Window;

  time(&t0);

  // Count the number of pixels in the object file to set g_bins. 
  objects = fopen(argv[1], "r");
  assert(objects!=NULL);
  printf("#Using %s as the input object file.\n", argv[1]);
  healpix = read_healpix_map(argv[1], &nside, coords, ordering); 
  NSIDE = count_Healpix_pixels(argv[1], healpix);
  //assert(NSIDE==nside);
  nside = NSIDE; //Joy and I were having issues with nside = 0 this is a work around

  // Count the number of lines in the bandpower file to set g_bands. 
  bandpowers = fopen(argv[2], "r");
  assert(bandpowers!=NULL);
  printf("#Using %s as the input bandpower file.\n", argv[2]);
  g_bands = object_count(bandpowers, kLinesPerObject);

  // Strip the directory structure from the filename if present
  if (strrchr(argv[2], '/')!=NULL) {
    strcpy(name, strrchr(argv[2], '/')+1);
  } else {
    strcpy(name, argv[2]);
  }

  sprintf(filename, "KL_%s", name);
  output_KL = fopen(filename, "w");
  sprintf(filename, "C_%s", name);
  output_C = fopen(filename, "w");
  sprintf(filename, "Fisher_%s", name);
  output_Fisher = fopen(filename, "w");
  sprintf(filename, "Window_%s", name);
  output_Window = fopen(filename, "w");
  
  printf("#Bins = %ld, Bands = %ld\n", g_bins, g_bands);
  fflush(NULL);

  C = (double *)malloc(g_bands*sizeof(double));
  C_change = (double *)malloc(g_bands*sizeof(double));
  C_start = (int *)malloc((g_bands)*sizeof(int));
  C_stop = (int *)malloc((g_bands)*sizeof(int));
  ra = (double *)malloc(g_bins*sizeof(double));
  dec = (double *)malloc(g_bins*sizeof(double));
  overdensity = (double *)malloc(g_bins*sizeof(double));
  data_covariance = (double *)malloc(g_bins*g_bins*sizeof(double));
  noise = (double *)malloc(g_bins*g_bins*sizeof(double));
  cos_angle = (double *)malloc(g_bins*g_bins*sizeof(double));
  signal = (double *)malloc(g_bands*g_bins*g_bins*sizeof(double));
  assert(C!=NULL);
  assert(C_change!=NULL);
  assert(C_start!=NULL);
  assert(C_stop!=NULL);
  assert(ra!=NULL);
  assert(dec!=NULL);
  assert(overdensity!=NULL);
  assert(data_covariance!=NULL);
  assert(noise!=NULL);
  assert(cos_angle!=NULL);
  assert(signal!=NULL);
  assert(fflush(NULL)==0);

  // Read the bandpower file to find C_start, C_stop, and initial value. 
  read_bandpower_file(bandpowers, C, C_start, C_stop);

  // Read the objects file to calculate the overdensities, solid angle, and pixel positions.  
  read_Healpix_file(overdensity, healpix, ra, dec, nside);

  // Calculate the covariance matrix based on the observed overdensities.
  time(&t2);
  calculate_Healpix_covariance(cos_angle, noise, data_covariance, ra, dec, overdensity);
  time(&t3);
  printf("#Calculated Healpix covariance.  Elapsed time = %lf seconds.\n", difftime(t3, t2));

  // Calculate the signal matrix using the angles between the pixels.
  time(&t2);
  signaltime = calculate_signal(signal, cos_angle, C_start, C_stop);
  time(&t3);
  printf("#Calculated signal.  Elapsed time = %lf seconds.\n", difftime(t3, t2));

  // KL-Compress
  time(&t2);
  KL_compression(overdensity, signal, noise, data_covariance, C, output_KL); 
  time(&t3);
  printf("#Calculated KL Compression.  Elapsed time = %lf seconds.\n", difftime(t3, t2));

  F = (double *)malloc(g_bands*g_bands*sizeof(double));
  A = (double *)malloc(g_bands*g_bins*g_bins*sizeof(double));
  B = (double *)malloc(g_bands*g_bins*g_bins*sizeof(double));
  average = (double *)malloc(g_bands*sizeof(double));
  model_covariance = (double *)malloc(g_bins*g_bins*sizeof(double));
  difference = (double *)malloc(g_bins*g_bins*sizeof(double));
  assert(F!=NULL);
  assert(A!=NULL);
  assert(B!=NULL);
  assert(average!=NULL);
  assert(model_covariance!=NULL);
  assert(difference!=NULL);

  // Begin iterative process to find C_l. 
  for (n=1; n<=kMaxIter; n++) {

    // Calculate the C_l and Fisher matrix from the signal and covariance matrices
    time(&t2);
    iterationtime += estimate_C(signal, model_covariance, data_covariance, noise, difference, average, A, B, F, C, C_start, C_stop, C_change, n, output_C, output_Fisher, output_Window);
    assert(fflush(NULL)==0);
    time(&t3);
    printf("#Calculated iteration %d.  Elapsed time = %lf seconds.\n", n, difftime(t3, t2));

  }

  time(&t1);

  // Print out time taken for the program and how much was done in the kernel. 
  printf("#End of program.  Elapsed time = %lf seconds.  Kernel is %lf percent of total time.\n", difftime(t1, t0), (signaltime+iterationtime)*100.0/(difftime(t1, t0)*1.0));
  assert(fflush(NULL)==0);

  return 0;
}
