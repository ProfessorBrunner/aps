#include "angular_power_spectrum.h"

int object_count(FILE *objects, int n)
     /* Count the number of objects in a file by counting the lines and dividing by the number of lines per object. 
	Ignore lines beginning with #. */
{
  int number_lines;
  char line[kMaxChars];
  
  number_lines = 0;
  while (fgets(line, sizeof(line), objects)!=NULL)
    if (line[0]!='#') number_lines++;
  rewind(objects);
  printf("#Done counting %d lines.\n", number_lines/n);
  
  return number_lines/n;
}

void tic(Timer *timer) {
  gettimeofday(&timer->start, NULL);
}

void toc(Timer *timer) {
  gettimeofday(&timer->stop, NULL);
  timer->elapsed = ((timer->stop.tv_sec - timer->start.tv_sec)*1000000 + 
    (timer->stop.tv_usec - timer->start.tv_usec) ) / 1000000.0;
}

void toc_print(Timer *timer) {
  toc(timer);
  printf("elapsed time: %f\n", timer->elapsed);
}

void save_raw_double_array(char *root, char *name, double *data, int len) {
  char filename[kMaxChars];
  int count;
  sprintf(filename, "%s/%s.dat", root, name);
  printf("#Writing test data %s to %s\n", name, filename);

  FILE *fp = fopen(filename, "wb");
  assert(fp);
  count = fwrite(data, sizeof(*data), len, fp);
  assert(count==len);
  fclose(fp);
}