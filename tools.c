#include "angular_power_spectrum.h"

int object_count(FILE *objects, int n)
     /* Count the number of objects in a file by counting the lines and dividing by the number of lines per object. 
	Ignore lines beginning with #. */
{
  int number_lines;
  char line[MAXCHARS];
  
  number_lines = 0;
  while (fgets(line, sizeof(line), objects)!=NULL)
    if (line[0]!='#') number_lines++;
  rewind(objects);
  printf("#Done counting %d lines.\n", number_lines/n);
  
  return number_lines/n;
}
