// SATID.cpp
#include <ctype.h>         // for isdigit in read_tle
#include "orbit_gray.h"    // for Date, Locate, Satellite, de2ra, etc.
using namespace std;

///////////// DECLARE GLOBAL VARIABLES ////////////////////////////////////////

double
   tle,      // epoch of elements in tle format
   ii,       // inclination, degrees
   om,       // right ascension of ascending node, degrees
   ec,       // eccentricity
   ww,       // argument of the perigee, degrees
   ma,       // mean anomaly, degrees
   nn,       // mean motion, revolutions/day
   uu,       // true longitude
   c2,       // bstar coefficient
   bstar;    // BSTAR drag term
char name[81];

//////////// DECLARE FUNCTIONS ////////////////////////////////////////////////

double sci(char *string);
void read_tle(char *file);
inline char *s_in(char *prompt, char *buffer);
void write_tle(char *file_out);

/////////////////// MAIN //////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  char id_file[20], in_file[20];
  char buf[81], buf1[81], buf2[81];
  double delt, delta_t, cpa, Perr, alpha, err1, err2;
  double v1[3], v2[3], vp1[3], vp2[3], delr[3];

  sprintf(in_file, "refer.tle");           // default reference file

  if(argc == 1)
    sprintf(id_file, "unid.txt");          // default unknown file
  if(argc == 2)
  {
    sprintf(id_file, argv[1]);
  }
  if(argc == 3)
  {
    sprintf(id_file, argv[1]);   // find this satellite (tle file)
    sprintf(in_file, argv[2]);   // in this archive of satellites
  }

  // read tle from id_file
  // orbit parameters are passed from read_tle as global variables
  // echo tle to screen
  read_tle(id_file);

  // id_sat comparison variables
  Date t1(tle);
  Satellite id_sat(t1.jd, ii, om, ec, ww, ma, nn, bstar);
  // flat projection of satellite direction unit vector, vp1
  cross(id_sat.rr, id_sat.vv, v1);
  cross(v1, id_sat.rr, v2);
  proj(v2, id_sat.vv, vp1);

  // set error bound
  err1 = 2;
  printf("   position error, degrees [%.1f", err1);
  if (strlen(s_in("]: ", buf))) err1 = atof(buf);
  err2 = 20;
  printf("track angle error, degrees [%.1f", err2);
  if (strlen(s_in("]: ", buf))) err2 = atof(buf);
  printf("\n");

  // read all TLE's in in_file, three line format, no echo
  FILE *fpi, *fpo;
  fpi = fopen(in_file, "r");
  // read in_sat
  while(fgets(name, 80, fpi))
  {
    if(strlen(name) == 70) continue;
    fgets(buf1, 80, fpi);        // line1
    fgets(buf2, 80, fpi);        // line2

    // first data line
    tle = atof(buf1 + 18);       // epoch in tle format at this point
    bstar = sci(buf1 + 53);
    ssn = atoi(buf1 + 2);        // satellite number
    sscanf(buf1 + 9, "%6s", &desig);      // international designation

    // second data line
    ii = atof(buf2 + 8);         // inclination, degrees
    om = atof(buf2 + 17);        // ascending node, degrees
      buf2[25] = '.';            // add decimal point to eccentricity
    ec = atof(buf2 + 25);        // eccentricity
    ww = atof(buf2 + 34);        // perigee, degrees
    ma = atof(buf2 + 43);        // mean anomaly, degrees
    // Make sure mean motion is null-terminated, since rev. no.
    // may immediately follow.
    buf[0] = buf2[63];
    buf2[63] = '\0';
    nn = atof(buf2 + 52);        // mean motion, revolutions/day
    buf2[63] = buf[0];

    // create sat
    Satellite sat(tle, ii, om, ec, ww, ma, nn, bstar);

    // advance sat vectors to id_sat epoch
    sat.delta_t(t1.jd);

    // advance elements to id epoch
    sat.rv2el(sat.rr, sat.vv);

    // satellite delta r vector, delr
    vmadd(id_sat.rr, sat.rr, delr, -1);

    // flat projection of delta r unit vector, delr
    cross(id_sat.rr, delr, v1);
    cross(v1, id_sat.rr, v2);
    proj(v2, delr, delr);

    // angle between delr and id_sat.vv, radians
    alpha = acos(dot(delr, id_sat.vv)/norm(id_sat.vv));

    // angle between position unit vectors, Perr, in radians.
    Perr = acos(dot(sat.rr, id_sat.rr)/(norm(sat.rr)*norm(id_sat.rr)));

    // magnitude of Perr in direction of id_sat.vv, radians
    delt = atan(tan(Perr) * cos(alpha));    // Napier's Rule
    // time of flight to cpa, seconds
    delta_t = 60*delt * norm(id_sat.rr) / norm(id_sat.vv);
    // closest point, radians
    cpa = asin(sin(alpha) * sin(Perr));     // Napier's Rule

    // flat projection of satellite direction unit vector, vp2
    cross(sat.rr, sat.vv, v1);
    cross(v1, sat.rr, v2);
    proj(v2, sat.vv, vp2);

    // angle between direction unit vectors, radians
    alpha = acos(dot(vp1, vp2));

    // deltas
    alpha  = acos(cos(alpha)/cos(delt));  // spherical law of cosines
    alpha /= de2ra;
    Perr  /= de2ra;

    // compare id_sat to in_sat using osculating elements (close enough)
    if(Perr < err1 && alpha < err2)
    {
       // if match is found

       // this is so we can write the updated tle to the DOS window
       tle = sat.tle;
       ii  = sat.xincl / de2ra;
       om  = sat.xnodeo / de2ra;
       ec  = sat.eo;
       ww  = sat.omegao / de2ra;
       ma  = sat.xmo / de2ra;
       nn  = sat.xno / nocon;
       c2  = sat.c2;

       // visually check match parameters using advanced mean elements
       printf("%s", name);
       write_tle("");
       printf("   position error  %4.1f\n", Perr);
       printf("track angle error  %4.1f\n\n", alpha);
       printf("       time error  %4.0f\n", delta_t);
       printf(" to closest point  %4.1f\n\n", cpa/de2ra);

       // output original tle to unid file
       fpo = fopen(id_file, "a");
       fprintf(fpo, "\n%s", name);
       fprintf(fpo, "%s", buf1);
       fprintf(fpo, "%s", buf2);
       fclose(fpo);
       s_in("\n[Next]", buf);
     }    // if match
  }     // while
  fclose(fpi);
  s_in("\n[Done]", buf);
  system(id_file);
}    // end main

////////////////// FUNCTIONS //////////////////////////////////////////////////

inline char *s_in(char *prompt, char *buffer)
{
  printf("%s", prompt);
  if (feof(stdin)) exit(0);
  gets(buffer);
  return buffer;
}

/* Converts quasi scientific notation to double.  Format:  SnSn where S is
either '-' or '\0' and n represents a sequence of 1 or more digits.  An
implied decimal point exists to the left of the left n.  Total length 27
chars max. Subroutine of read_tle*/
double sci(char *string)
{
   char buf[30], *bufptr;

   bufptr = buf;

   if (string[1] == '\0')
      return 0.;

   /* get significand */
   if (*string == '-')
      *bufptr++ = '-';
   *bufptr = '.';
   while (isdigit(*++bufptr = *++string))
      ;   /* copy significand */
   *bufptr++ = 'E';

   /* get exponent */
   if (*string == '\0')   /* no sign on exponent */
      ++string;

   strcpy(bufptr, string);      /* copy exponent */
   return atof(buf);
}

// smart tle reader, reads first tle in file, no conversion of elements
void read_tle(char *file)
{
  FILE *fp;
  fp = fopen(file, "r");
  char buf1[81], buf2[81];
  buf1[0] = '0';
  do
  {
    if(*buf1 == '1' && strlen(buf1) == 70)
    {
      fgets(buf2, 80, fp);
      if(*buf2 == '2' && strlen(buf2) == 70)
      {
        fclose(fp);

        // first data line, epoch and bstar
        printf("%s", buf1);         // print first line on screen
        tle = atof(buf1 + 18);      // epoch in tle format at this point
        bstar = sci(buf1 + 53);

        // second data line
        printf("%s\n", buf2);        // print second line on screen
        ii = atof(buf2 + 8);         // inclination, degrees
        om = atof(buf2 + 17);        // ascending node, degrees
          buf2[25] = '.';            // add decimal point to eccentricity
        ec = atof(buf2 + 25);        // eccentricity
        ww = atof(buf2 + 34);        // perigee, degrees
        ma = atof(buf2 + 43);        // mean anomaly, degrees
        // Make sure mean motion is null-terminated, since rev. no.
        // may immediately follow.
        buf2[63] = '\0';
        nn = atof(buf2 + 52);        // mean motion, revolutions/day
        return;
      }
    }
  }while(fgets(buf1, 80, fp));

  printf("\nNo TLE found in file %s\n", file);
  s_in("[exit]", buf1);
  exit(0);
}

// write TLE to output file and to screen
void write_tle(char *file_out)
{
  char ec_string[9];
  char bstar_string[13];
  char bstar_fract[7];
  char bstar_exp[3];
  char line1[81], line2[81];

  sprintf(ec_string, "%.7f", ec);
  ec_string[0] = ec_string[2];
  ec_string[1] = ec_string[3];
  ec_string[2] = ec_string[4];
  ec_string[3] = ec_string[5];
  ec_string[4] = ec_string[6];
  ec_string[5] = ec_string[7];
  ec_string[6] = ec_string[8];

  sprintf(bstar_string, "%12.4e", bstar*10);
  bstar_fract[0] = bstar_string[0]; // sign
  bstar_fract[1] = bstar_string[1];
  bstar_fract[2] = bstar_string[3];
  bstar_fract[3] = bstar_string[4];
  bstar_fract[4] = bstar_string[5];
  bstar_fract[5] = bstar_string[6];
  bstar_fract[6] = '\0';
  bstar_exp[0] = bstar_string[8];
  bstar_exp[1] = bstar_string[11];
  bstar_exp[2] = '\0';

  double xns = 2160 * bstar * nn * c2;

  sprintf(line1, "1 %05dU %-8s %014.8f %.8f  00000-0 %6s%2s 0    00"
           ,ssn, desig, tle, xns, bstar_fract, bstar_exp);
  ccksum(line1);
  sprintf(line2, "2 %05d %8.4lf %8.4lf %.7s %8.4lf %8.4lf %11.8lf    00",
           ssn, ii, om, ec_string, ww, ma, nn);
  ccksum(line2);

  printf("%s\n", line1);
  printf("%s\n", line2);
}


