// FIX96.cpp
// by Scott Campbell campbel7@hughes.net
//

/******************************************************************************
*                                                                             *
*                                                                             *
*                                                                             *
*******************************************************************************/

#include <cstdlib>
#include <cstdio>
#include <string.h>
using namespace std;

void ccksum(char *line);

/////////////////// MAIN ///////////////////////////////////////////////////////
main(void)
{
  char file[81], line1[81], line2[81];
  char buf[81], desig[6];
  double tle;
  int ssn;
  FILE *fpi, *fpo;

  system("copy int96tles.tle list.txt");

  // open record file
  fpi = fopen("list.txt", "r");
  // create batch file as output
  fpo = fopen("int96tles.tle", "w");

  // read all the lines
  while(fgets(line1, 80, fpi))
  {
    if(*line1 == '1' && strlen(line1) == 70)
    {
      fgets(line2, 80, fpi);
      if(*line2 == '2' && strlen(line2) == 70)
      {
        // first data line
        ssn = atoi(line1 + 2);
        if((ssn - 96000) > 0)
        {
          sprintf(desig, "96%3dA", ssn - 95500);
        }
        else
        {
          sscanf(line1 + 9, "%6s", &desig);
        }
        tle = atof(line1 + 18);
      }
    }
    sprintf(line1, "1 %05dU %-8s %014.8f 0.00000000  00000-0  00000-0 0    00"
                      ,ssn, desig, tle);
    ccksum(line1);
    fprintf(fpo, "%s\n", line1);
    fprintf(fpo, "%s", line2);



  } // end of all files in directory

  fclose(fpi);
  fclose(fpo);
  system("del list.txt");
}    // end main

void ccksum(char *line)
{
    int cksum, i;

    cksum = 0;
    for (i = 0; i < 68; i++) {
        char ich;
        ich = line[i];
        if (ich > '0' && ich <= '9') cksum += ich - '0';
        if (ich == '-') cksum += 1;
    }
    line[68] = (cksum % 10) + '0';
    line[69] = '\0';
}

