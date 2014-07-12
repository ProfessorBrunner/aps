#include "angular_power_spectrum.h"
//#include "new_angular_power_spectrum.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

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
  Timer time_total;
  Timer time_find_c_iteration;
  Timer time_function;
  tic(&time_total);

  int i, n, *C_start, *C_stop;
  long nside, NSIDE;
  float *healpix;
  double *overdensity, *data_covariance, *noise, *cos_angle, *difference, *A, *B, *C, *C_change; 
  double *ra, *dec, *F, *average, *model_covariance, *signal, signaltime, iterationtime = 0;
  char ordering[10], coords[1], filename[kMaxChars], output_root[kMaxChars], 
       name[kMaxChars], test_root[kMaxChars], *char_position;
  FILE *objects, *bandpowers, *output_KL, *output_C, *output_Fisher, *output_Window;

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
  assert(strlen(argv[2])<kMaxChars);
  strcpy(output_root, argv[2]);
  char_position = strrchr(output_root, '/');
  if (char_position) {
    strcpy(name, char_position+1);
    *char_position = '\0';
  } else {
    strcpy(name, output_root);
    output_root[0] = '.';
    output_root[1] = '\0';
  }


  struct stat st = {0};

  sprintf(filename, "%s/output", output_root);
  if (stat(filename, &st) == -1) mkdir(filename, 0766);

# ifdef APS_OUTPUT_TEST
  if (argc > 3) {
    sprintf(test_root, "%s/test_shared_%s_%s", output_root, name, argv[3]);
  } else {
    sprintf(test_root, "%s/test_shared_%s", output_root, name);
  }
  char_position = strrchr(test_root, '.');
  if (char_position) *char_position = '\0';
  if (stat(test_root, &st) == -1) mkdir(test_root, 0766);
# endif

  sprintf(filename, "%s/output/KL_%s", output_root, name);
  printf("%s\n", filename);
  output_KL = fopen(filename, "w");
  assert(output_KL);

  sprintf(filename, "%s/output/C_%s", output_root, name);
  output_C = fopen(filename, "w");
  assert(output_C);

  sprintf(filename, "%s/output/Fisher_%s", output_root, name);
  output_Fisher = fopen(filename, "w");
  assert(output_Fisher);

  sprintf(filename, "%s/output/Window_%s", output_root, name);
  output_Window = fopen(filename, "w");
  assert(output_Window);


  
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
  tic(&time_function);
  calculate_Healpix_covariance(cos_angle, noise, data_covariance, ra, dec, overdensity);
  toc(&time_function);
  printf("#Calculated Healpix covariance.  Elapsed time = %g seconds.\n", time_function.elapsed);

  // Calculate the signal matrix using the angles between the pixels.
  tic(&time_function);
  signaltime = calculate_signal(signal, cos_angle, C_start, C_stop);
  toc(&time_function);
  printf("#Calculated signal.  Elapsed time = %g seconds.\n", time_function.elapsed);
# ifdef APS_OUTPUT_TEST
  save_raw_double_array(test_root, "overdensity", overdensity, g_bins);
  for (i = 0; i < g_bands; ++i){
    sprintf(filename, "signal%d", i);
    save_raw_double_array(test_root, filename, (signal + g_bins*g_bins*i), g_bins*g_bins);
  }
# endif

  // KL-Compress
  tic(&time_function);
  KL_compression(overdensity, signal, noise, data_covariance, C, output_KL, test_root); 
  toc(&time_function);
  printf("#Calculated KL Compression.  Elapsed time = %g seconds.\n", time_function.elapsed);
# ifdef APS_OUTPUT_TEST
  save_raw_double_array(test_root, "kl_overdensity", overdensity, g_bins);
  for (i = 0; i < g_bands; ++i){
    sprintf(filename, "kl_signal%d", i);
    save_raw_double_array(test_root, filename, (signal + g_bins*g_bins*i), g_bins*g_bins);
  }
  save_raw_double_array(test_root, "kl_noise", signal, g_bins*g_bins);
  save_raw_double_array(test_root, "kl_covariance_data", signal, g_bins*g_bins);
# endif

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
    tic(&time_find_c_iteration);
    iterationtime += estimate_C(signal, model_covariance, data_covariance, noise, difference, average, A, B, F, C, C_start, C_stop, C_change, n, output_C, output_Fisher, output_Window);
    assert(fflush(NULL)==0);
    toc(&time_find_c_iteration);
    printf("#Calculated iteration %d.  Elapsed time = %g seconds.\n", n, time_find_c_iteration.elapsed);
#   ifdef APS_OUTPUT_TEST
    //difference doesn't change
    if (n == 1) save_raw_double_array(test_root, "difference", difference, g_bins*g_bins);
    sprintf(filename, "iter_%d_covariance_model", n);
    save_raw_double_array(test_root, filename, model_covariance, g_bins*g_bins);
    sprintf(filename, "iter_%d_average", n);
    save_raw_double_array(test_root, filename, average, g_bands);
    sprintf(filename, "iter_%d_A", n);
    save_raw_double_array(test_root, filename, A, g_bins*g_bins*g_bands);
    sprintf(filename, "iter_%d_B", n);
    save_raw_double_array(test_root, filename, B, g_bins*g_bins*g_bands);
    sprintf(filename, "iter_%d_F", n);
    save_raw_double_array(test_root, filename, F, g_bands*g_bands);
    sprintf(filename, "iter_%d_C", n);
    save_raw_double_array(test_root, filename, C, g_bands);

#   endif
  }

  toc(&time_total);
  // Print out time taken for the program and how much was done in the kernel. 
  //printf("#End of program.  Elapsed time = %g seconds.  Kernel is %lf percent of total time.\n", difftime(t1, t0), (signaltime+iterationtime)*100.0/(difftime(t1, t0)*1.0));
  printf("#End of program.  Elapsed time = %g seconds.\n", time_total.elapsed);
  assert(fflush(NULL)==0);


  return 0;
}
