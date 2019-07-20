// RV2TEL_KM.cpp
// by Scott Campbell campbel.7@hotmail.com

/******************************************************************************
*  Supply a file named vec.txt with the position vector in units of km        *
*  and the velocity vector in units of km / sec.  Run the                     *
*  program and follow prompts.  Press ENTER to accept default epoch values.   *
*******************************************************************************/

#include <ctype.h>     // for isdigit in cksum
#include "orbit_gray.h"
using namespace std;


/*####################### UNITS & CONVENTIONS #########################

Unless otherwise indicated, throughout this program
quantities are measured in the following units:

angles        radians
length        equatorial earth radii (1 unit = xkmper km)
velocity      equatorial earth radii / minute
South latitudes are negative.
East longitude is positive, west is negative.

**********************************************************************/


///////////// declare global variables ////////////////////////////////////////

double tle, ii, om, ec, ww, ma, nn;
char buf[21], name[20],
     sn[] = "12345", dsg[] = "12345A", bstar[] = "10000-4";
double rr[3], vv[3];

//////////// declare functions ////////////////////////////////////////////////

void get_vec(char *file_in);
int cksum(char *line);
void write_tle(char *file_out);
char *s_in(char *prompt, char *buffer);

/////////////////// MAIN //////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  char file_in[20], file_out[20];


  if(argc == 2)
  {
    sprintf(file_in, argv[1]);
    sprintf(file_out, argv[1]);
  }
  else if(argc == 3)
  {
    sprintf(file_in, argv[1]);
    sprintf(file_out, argv[2]);
  }
  else          // default input/output files
  {
    sprintf(file_in, "VEC_KM.txt");
    sprintf(file_out, "VEC_KM.txt");
  }

  s_in("name           :", name);

  printf("ssn     [%5s", sn);
  if (strlen(s_in("]: ", buf))) sprintf(sn, "%5s", buf);

  printf("desig  [%6s", dsg);
  if (strlen(s_in("]: ", buf))) sprintf(dsg, "%6s", buf);

  printf("bstar [%7s", bstar);
  if (strlen(s_in("]: ", buf))) sprintf(bstar, "%7s", buf);

  get_vec(file_in);
  printf("Enter TLE Epoch\n");
  Date t1;
  t1.input();
  tle =  t1.tle;
  Satellite sat(tle, rr, vv, 0.);

  // satellite mean elements
  ii = sat.xincl / de2ra,     // inclination
  om = sat.xnodeo / de2ra,    // right ascension of ascending node
  ec = sat.eo,
  ww = sat.omegao / de2ra,    // argument of the perigee
  ma = sat.xmo / de2ra,       // mean anomaly
  nn = sat.xno / nocon;       // mean motion

  write_tle(file_out);
  s_in("\n[exit]", buf);
  return 0;
}

////////////////// end main ///////////////////////////////////////////////////

inline char *s_in(char *prompt, char *buffer)
{
  printf("%s", prompt);
  if (feof(stdin)) exit(0);
  gets(buffer);
  return buffer;
}

// load state vectors from input file_in
void get_vec(char *file_in)
{
  FILE *fp;
  if ((fp = fopen(file_in, "r")) == NULL)
  {
    printf("can't open %s\n", file_in);
    s_in("[exit]", buf);
    exit(0);
  }

  for(int i = 0; i < 6; i++)
  {
    fgets(buf, 20, fp);
    if(i < 3)  rr[i] = atof(buf) / 6378.135;
    if(i >= 3)  vv[i - 3] = atof(buf) * 60 / 6378.135;
  }
  fclose(fp);
}

// write TLE to output file and to screen
void write_tle(char *file_out)
{
  char ec_string[9];
  char line1[70];
  char line2[70];

  sprintf(ec_string, "%.7f", ec);
  ec_string[0] = ec_string[2];
  ec_string[1] = ec_string[3];
  ec_string[2] = ec_string[4];
  ec_string[3] = ec_string[5];
  ec_string[4] = ec_string[6];
  ec_string[5] = ec_string[7];
  ec_string[6] = ec_string[8];

  double xndt2o = 0;

  sprintf(line1, "1 %-.5sU %-.6s   0%13.8f %.8f  00000-0  %7s 0    00"
           ,sn, dsg, tle, xndt2o, bstar);
   // change tle format character from 0%13.8f to %14.8f on 1/1/2010
  ccksum(line1);
  sprintf(line2, "2 %-.5s %8.4lf %8.4lf %.7s %8.4lf %8.4lf %11.8lf    00",
    sn, ii, om, ec_string, ww, ma, nn);
  ccksum(line2);
  FILE *fp;
  fp = fopen(file_out, "a");
  fprintf(fp, "\n%s\n", name);
  fprintf(fp, "%s\n", line1);
  fprintf(fp, "%s\n", line2);
  fclose(fp);
  printf("\n%s\n", name);
  printf("%s\n", line1);
  printf("%s\n", line2);
}

