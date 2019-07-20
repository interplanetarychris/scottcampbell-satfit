// SST.CPP
#include <ctype.h>     // for isdigit in read_tle
#include "sst.h"       // for Date, Locate, Satellite

using namespace std;

///////////// DECLARE GLOBAL VARIABLES /////////////////////////////////////////

double
   jd,       // epoch of elements in julian date
   tle,      // epoch of elements in tle format
   ii,       // inclination, degrees
   om,       // right ascension of ascending node, degrees
   ec,       // eccentricity
   ww,       // argument of the perigee, degrees
   uu,       // longitude
   ma,       // mean anomaly, degrees
   nn,       // mean motion, revolutions/day
   xns,      // ndot
   c2,       // internal drag term
   bstar;    // BSTAR drag term


//////////// DECLARE FUNCTIONS /////////////////////////////////////////////////

double sci(char *string);
void write_tle(FILE *file);

/////////////////// MAIN ///////////////////////////////////////////////////////

char buf[80];
char name[80], line1[80], line2[80];
FILE *output;

int main(int argc, char *argv[])
{
   FILE *input;

   int year_limit, month_limit, day_limit;

   /* Get the input/output file name and open it */
   if(argc == 1)
   {
     input = fopen("advance.txt","r");
     output = fopen("advance.out","a");
   }
   if(argc == 2)
   {
     input = fopen(argv[1],"r");
     output = fopen("advance.out","a");
   }
   if(argc == 3)
   {
     input = fopen(argv[1],"r");
     output = fopen(argv[2],"a");
   }

   fgets(buf, 80, input); // only one time limit line
   sscanf(buf, "%d %d %d", &year_limit, &month_limit, &day_limit);
   Date t1(year_limit, month_limit, day_limit, 0., 0., 0.);

   while(fgets(buf, 80, input))     // line 1 or name
   {
      name[0] = '\0';
      if (buf[0] != '1')            // it's name
      {
         strcpy(name, buf);
         if (fgets(buf, 80, input) == NULL)  // line 1
            goto done;
      }
      strcpy(line1, buf);
      // first data line
      ssn = atoi(line1 + 2);        // satellite number
      sscanf(line1 + 9, "%6s", &desig);      // international designation
      tle = atof(line1 + 18);       // epoch in tle format at this point
      xns = atof(line1 + 33);
      bstar = sci(line1 + 53);
      if (fgets(buf, 80, input) == NULL)     // line 2
         goto done;
      strcpy(line2, buf);
      // second data line
      ii = atof(line2 + 8);         // inclination, degrees
      om = atof(line2 + 17);        // ascending node, degrees
        line2[25] = '.';            // add decimal point to eccentricity
      ec = atof(line2 + 25);        // eccentricity
      ww = atof(line2 + 34);        // perigee, degrees
      ww = mod(ww, 360);
      ma = atof(line2 + 43);        // mean anomaly, degrees
      // Make sure mean motion is null-terminated, since rev. no.
      // may immediately follow.
      buf[0] = line2[63];
      line2[63] = '\0';
      nn = atof(line2 + 52);        // mean motion, revolutions/day
      line2[63] = buf[0];           // restore line2[63]

      printf("\n");
      Satellite sat(tle, ii, om, ec, ww, ma, nn, bstar);
      printf("%s", name);
      sat.print_el();

      // SST
      sat.delta_sst(t1.jd, 1);      // 1 = advance to node
      sat.rv2el(sat.rr, sat.vv);

      // write_tle variables
      tle = sat.tle;
      ii = sat.xincl/de2ra;
      om = sat.xnodeo/de2ra;
      ec = sat.eo;
      ww = sat.omegao/de2ra;
      ma = sat.xmo/de2ra;
      c2 = sat.c2;
      nn = sat.xno/nocon;
      bstar = sat.bstar;

      // advance output to file and screen
      printf("\nadvanced to %s", name);
      fprintf(output, "%s", name);
      write_tle(output);

   }  // end input while

done:

   fclose(input);
   fclose(output);
   s_in("\n[done]", buf);
   return(0);
}


////////////////// FUNCTIONS //////////////////////////////////////////////////

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

// write TLE to output file and to screen
void write_tle(FILE *file_out)
{
  char ec_string[9];
  char xns_string[12];
  char bstar_string[13];
  char bstar_fract[7];
  char bstar_exp[3];

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

  xns = 2160 * bstar * nn * c2;

  sprintf(xns_string, "%.8f", xns);
  if(xns_string[0] == '-')
  {
    xns_string[1] = xns_string[2];
    xns_string[2] = xns_string[3];
    xns_string[3] = xns_string[4];
    xns_string[4] = xns_string[5];
    xns_string[5] = xns_string[6];
    xns_string[6] = xns_string[7];
    xns_string[7] = xns_string[8];
    xns_string[8] = xns_string[9];
    xns_string[9] = xns_string[10];
    xns_string[10] = xns_string[11];
  }

  sprintf(line1, "1 %05dU %-8s %014.8f %10s  00000-0 %6s%2s 0    00"
           ,ssn, desig, tle, xns_string, bstar_fract, bstar_exp);
  ccksum(line1);
  sprintf(line2, "2 %05d %8.4lf %8.4lf %.7s %8.4lf %8.4lf %11.8lf    00",
           ssn, ii, om, ec_string, ww, ma, nn);
  ccksum(line2);


  fprintf(file_out, "%s\n", line1);
  fprintf(file_out, "%s\n", line2);

  printf("%s\n", line1);
  printf("%s\n", line2);
}

