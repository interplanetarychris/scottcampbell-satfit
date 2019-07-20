// TrendA.cpp
// by Scott Campbell campbel7@hughes.net
//

/******************************************************************************
*                                                                             *
*                                                                             *
*                                                                             *
*******************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "date.h"   // for time object
using namespace std;


inline char *s_in(char *prompt, char *buffer)
{
  printf("%s", prompt);
  if (feof(stdin)) exit(0);
  gets(buffer);
  return buffer;
}

/////////////////// MAIN ///////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  char line1[81], line2[81], file_in[81], file_out[81];
  char element[2], buf[81];
  char *str;
  double tle, zero, ii, om, ec, ww, ma, nn;
  int i, ssn;
  FILE *fp, *fpo, *fps;

  ssn = 27711;

  top:

  if(argc == 2)
  {
    str = strstr(argv[1], ".");
    sscanf(str-5, "%s", &buf);
    ssn = atoi(buf);
    sprintf(file_in, "history/%05d.txt", ssn);
    printf("\nopening %05d.txt\n", ssn);
  }
  else
  {
    printf("\nsatellite ssn [%d", ssn);
    if (strlen(s_in("]: ", buf)))
      ssn = atoi(buf);
    if ((*buf & 0x5f) == 'Q') exit(0);
    sprintf(file_in, "history/%05d.txt", ssn);
  }
  argc = 1;

  i = 0;
  if(fp = fopen(file_in, "r"))
  {

    again:
    printf("\n(A)ll (I)ncl (O)mega (E)cc (P)erigee (M)otion (H)istory (Q)uit: ");
    s_in("", element);
    if ((*element & 0x5f) == 'Q') exit(0);
    if ((*element & 0x5f) == 'A') sprintf(file_out, "trendA.csv");
    else sprintf(file_out, "trend.csv");
    fpo = fopen(file_out, "w");
    if ((*element & 0x5f) == 'H')
    {
      fps = fopen("batch.inp", "w");
      fprintf(fps, "%s", "y");
      fclose(fps);
      sprintf(buf, "satfit %s ba<batch.inp", file_in);
      system(buf);
      sprintf(buf, "DEL batch.inp");
      system(buf);
      goto again;
    }
    fprintf(fpo, "time");
    if ((*element & 0x5f) == 'I' || (*element & 0x5f) == 'A')
          fprintf(fpo, ",inc");
    if ((*element & 0x5f) == 'O' || (*element & 0x5f) == 'A')
          fprintf(fpo, ",omega");
    if ((*element & 0x5f) == 'E' || (*element & 0x5f) == 'A')
          fprintf(fpo, ",ecc");
    if ((*element & 0x5f) == 'P' || (*element & 0x5f) == 'A')
          fprintf(fpo, ",perigee");
    if ((*element & 0x5f) == 'M' || (*element & 0x5f) == 'A')
          fprintf(fpo, ",motion");
    fprintf(fpo, "\n");

    while(fgets(line1, 80, fp))
    {
      if(*line1 == '1' && strlen(line1) == 70)
      {
        fgets(line2, 80, fp);
        if(*line2 == '2' && strlen(line2) == 70)
        {
          // first data line
          tle = atof(line1 + 18);      // epoch in tle format at this point
          Date t1(tle);
          tle = t1.jd;
          if(i == 0)
          {
            zero = tle;
            i++;
          }
          tle -= zero;

          // second data line
          ii = atof(line2 + 8);         // inclination, degrees
          om = atof(line2 + 17);        // ascending node, degrees
            line2[25] = '.';            // add decimal point to eccentricity
          ec = atof(line2 + 25);        // eccentricity
          ww = atof(line2 + 34);        // perigee, degrees
          ma = atof(line2 + 43);        // mean anomaly, degrees
          // Make sure mean motion is null-terminated, since rev. no.
          // may immediately follow.
          line2[63] = '\0';
          nn = atof(line2 + 52);        // mean motion, revolutions/day
        }

        fprintf(fpo, "%f", tle);
        if ((*element & 0x5f) == 'I' || (*element & 0x5f) == 'A')
              fprintf(fpo, ",%f", ii);
        if ((*element & 0x5f) == 'O' || (*element & 0x5f) == 'A')
              fprintf(fpo, ",%f", om);
        if ((*element & 0x5f) == 'E' || (*element & 0x5f) == 'A')
              fprintf(fpo, ",%f", ec);
        if ((*element & 0x5f) == 'P' || (*element & 0x5f) == 'A')
              fprintf(fpo, ",%f", ww);
        if ((*element & 0x5f) == 'M' || (*element & 0x5f) == 'A')
              fprintf(fpo, ",%f", nn);
        fprintf(fpo, "\n");

      }
    }
    fclose(fpo);
    system(file_out);
    sprintf(buf, "del %s", file_out);
    system(buf);
  }
  else
  {
    printf("no file named %s\n", file_in);
  }
  fclose(fp);

  goto top;

}    // end main


