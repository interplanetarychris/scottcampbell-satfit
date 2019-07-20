/* ELCORD.CPP
by Paul S. Hirose 1990 Sep 30
Main file for orbital element correction program.
*/

/* master header */
#define MAIN
#include <string.h>
#include "elcor.h"

/*########################### DATA ###########################*/

/* calendar - days in full months */
static int calndr[13] = {0, 31, 59, 90, 120, 151, 181, 212, 243,
   273, 304, 334, 365};

static int batch = 0;

/*##################### LOCAL FUNCTIONS #####################*/

static void bye();      /* free heap, exit */
static void cpy0p(char *dest, char *src, int n);   /* string copy with '0' padding */
static int leapyr(int y);   /* returns true if year is leap */
static char *stoup();   /* convert string to upper case */
static void write_elem(char *file); /* write out elements to file */

/*########################### CODE ###########################*/

int main(int argc, char **argv)
{
   printf("\nELCORD\n");
   int size;
   struct elemen *elemp;
   unsigned u;
   char file[81], file_in[81];              // input file strings declared
   char tle1[81], tle2[81];
   char *cp;
   int iflag;
   FILE *fp;
   out = 0;
   newrms = 1000.;

   if(argc == 1)
     sprintf(file, "sat.txt");            //  default input/output file

   if(argc == 2)                            // command line file or default?
     sprintf(file, argv[1]);                // command line input/output file

   if (argc == 3 && strcmp(argv[2],"ba") == 0)
   {
     batch = 1;
     sprintf(file, argv[1]);
   }
   if(strchr(file, '\\'))
     sprintf(file, "%s", strrchr(file, '\\')+1);

   /* Load master obs sites file */
   getsites();

   /* Load input data */
   getobs(file);
   write_elem("");
   sprintf(tle1, "%s", line1);  // copy original line1 to tle1
   sprintf(tle2, "%s", line2);  // copy original line2 to tle2

   goto accept;

select:
   /* list the names of the orbital elements, ask which to refine */
   elemp = elarr;
   for (u = 0; u < NUMEL; u++, elemp++)
      printf("%u %s\n", u, elemp->elname);
   printf("Which ones to refine (up to %u)", min_val(nobs * 2, NUMEL));
   s_in(":", buf);
   tok();
   p = 0;
   while (1) {
      cp = *tokp;
      if (cp == NULL) break;
      tomod[p++] = atoi(cp);
      tokp++;
   }

   if ((PX2P = (double *)malloc(sizeof(double) * 2 * p * p)) == NULL)
          { printf("malloc failure\n"); exit(2); }

   if ((ATW = (double *)malloc(sizeof(double) * p * (n + 1))) == NULL)
          { printf("malloc failure\n"); exit(2); }

   if ((AT = (double *)malloc(sizeof(double) * n)) == NULL)
          { printf("malloc failure\n"); exit(2); }

diffcorr:
   printf("\nold rms %.5f - ", oldrms);
   do
   {
      for(u = 0; u < p; ++u)
          deriv(u);
      oldrms = newrms;
      lsqr();
      prnres();   /* compute residuals */
   }while(oldrms > newrms);

   printf("new rms %.5f\n", newrms);
   out = 0;
   write_elem("");

   /* begin the main loop of the program */

accept:
   /* Accept a new command */
   s_in("\nEnter command (S)el  (D)iff  (C)hg  (R)es  (B)atch  (W)rt, (Q) : ", buf);

   if ((*buf & 0x5f) == 'S') goto select;
   if ((*buf & 0x5f) == 'D') goto diffcorr;
   if ((*buf & 0x5f) == 'C') goto change;
   if ((*buf & 0x5f) == 'B') goto batch;
   if ((*buf & 0x5f) == 'I') goto init;       // hidden
   if ((*buf & 0x5f) == 'R')
   {
      out = 1;
      goto residuals;
   }
   if ((*buf & 0x5f) == 'E') goto edit_file;  // hidden
   if ((*buf & 0x5f) == 'W') goto write_el;
   if ((*buf & 0x5f) == 'Q') bye();
   goto accept;

change:
   /* print orbital elements and accept changes */
   iflag = 0;
   elemp = elarr;
   for (u = 0; u < NUMEL; u++, elemp++)
      el_in(u, elemp, elemp->elp, &iflag);

   printf("mean anomaly %.4f\n", (qmo - omegao) * ra2de);
   if (iflag)      /* you modified at least one element */
      goto residuals;
   goto accept;

residuals:
   printf("\nResiduals:\n");
   prnres();   /* compute residuals */
   out = 0;
   write_elem("");
   goto accept;

edit_file:
   system(file);
   goto accept;

batch:
   if(((1-eo)*pow(xke/ xno, 2./3)-1)*6378.135 < 2000)
   {
     // create batch file
     fp = fopen("batch.inp", "w");
     fprintf(fp, "%s\n", "s");
     fprintf(fp, "%d %d\n", 0, 4);
     fprintf(fp, "%s\n", "d");
     fprintf(fp, "%s\n", "s");
     fprintf(fp, "%d %d\n", 3, 4);
     fprintf(fp, "%s\n", "d");
     fprintf(fp, "%s\n", "s");
     fprintf(fp, "%d %d %d\n", 0, 3, 4);
     fprintf(fp, "%s\n", "d");
     fprintf(fp, "%s\n", "s");
     fprintf(fp, "%d %d\n", 5, 6);
     fprintf(fp, "%s\n", "d");
     fprintf(fp, "%s\n", "s");
     fprintf(fp, "%d %d %d %d\n", 0, 3, 4, 5);
     fprintf(fp, "%s\n", "d");
     fprintf(fp, "%s\n", "d");
     fprintf(fp, "%s\n", "s");
     fprintf(fp, "%d %d %d %d %d\n", 0, 3, 4, 5, 6);
     fprintf(fp, "%s\n", "d");
     fprintf(fp, "%s\n", "d");
     fprintf(fp, "%s\n", "s");
     fprintf(fp, "%d %d\n", 0, 1);
     fprintf(fp, "%s\n", "d");
     fprintf(fp, "%s\n", "s");
     fprintf(fp, "%d\n", 2);
     fprintf(fp, "%s\n", "d");
     fprintf(fp, "%s\n", "w");
     fprintf(fp, "%s\n", "u");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "q");
     fclose(fp);
     sprintf(file_in, "elcord %s ba<batch.inp", file);
     system(file_in);
     sprintf(file_in, "DEL batch.inp");
     system(file_in);
     sprintf(file_in, "elcord %s", file);
     system(file_in);
     exit(0);
   }
   else
   {
     goto init;
   }


init:
   // create batch file
   fp = fopen("batch.inp", "w");
   fprintf(fp, "%s\n", "s");
   fprintf(fp, "%d %d\n", 0, 4);
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "s");
   fprintf(fp, "%d %d\n", 0, 1);
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "s");
   fprintf(fp, "%d %d %d\n", 0, 1, 3);
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "s");
   fprintf(fp, "%d %d\n", 5, 6);
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "s");
   fprintf(fp, "%d %d %d %d\n", 0, 1, 3, 4);
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "s");
   fprintf(fp, "%d %d %d %d %d\n", 0, 1, 3, 4, 5);
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "s");
   fprintf(fp, "%d %d %d %d %d %d\n", 0, 1, 3, 4, 5, 6);
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "d");
   fprintf(fp, "%s\n", "w");
   fprintf(fp, "%s\n", "u");
   fprintf(fp, "%s\n", "q");
   fprintf(fp, "%s\n", "q");
   fclose(fp);
   sprintf(file_in, "elcord %s ba<batch.inp", file);
   system(file_in);
   sprintf(file_in, "DEL batch.inp");
   system(file_in);
   sprintf(file_in, "elcord %s", file);
   system(file_in);
   exit(0);

write_el:

   printf("\n(U)pdate (V)iew (A)ppend (H)istory (O)riginal (Q)uit");
   s_in(": ", buf);

   if ((*buf & 0x5f) == 'V') write_elem("");
   if ((*buf & 0x5f) == 'A') write_elem(file);
   if ((*buf & 0x5f) == 'O')
   {
     printf("\n%s\n", tle1);
     printf("%s\n", tle2);
   }
   if ((*buf & 0x5f) == 'H')
   {
      FILE *fp;
      sprintf(file_in, "history/%5d.txt", sat_ncat);
      if(sat_ncat < 10000) file_in[8] = '0';
      if ((fp = fopen(file_in, "r")) == NULL)
      {
        printf("NO HISTORY FOR THIS SAT");
        goto accept;
      }
      while(fgets(buf, 80, fp))   // get a line from ssn.txt
        printf("%s", buf);
      fclose(fp);
   }
   if ((*buf & 0x5f) == 'U')
   {
      char iod_line[100][81];
      FILE *fpi, *fpo;
      fpi = fopen(file, "r");

      nobs = 0;

      // find observation entries and write to iod_line matrix
      while(fgets(iod_line[nobs], 80, fpi))
      if(strlen(iod_line[nobs]) > 58 && strlen(iod_line[nobs]) != 70) nobs++;

      fclose(fpi);

      out = 1;
      write_elem(file);

      fpo = fopen(file, "a");
      fprintf(fpo, "\n");
      for(int i = 0; i < nobs; i++)
      {
         fprintf(fpo, "%s", iod_line[i]);
      }
      out = 0;
      fclose(fpo);
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto write_el;
}

//////////////// end main /////////////////////////////////

static double dint(double x)
{
   if (x >= 0.0)
      return floor(x);
   else
      return ceil(x);
}

/* Returns value (in minutes) of string of form "hhmm:ss.s...". */
double atomin(char *string)
{
   double time;
   char *ptr, c, buf[5];

   ptr = string;

   /* point to the colon (or terminating null if no colon) */
   while ((c = *ptr) != '\0' && c != ':')
      ptr++;

   /* process seconds */
   if (c) {      /* string had a colon */
      time = atof(ptr + 1) / 60.;
      *ptr = '\0';      /* turn ':' to '\0' */
   } else
      time = 0.;

   cpy0p(buf, string, 5);      /* buf[] = "hhmm" */

   time += atof(buf + 2);      /* add minutes */
   buf[2] = '\0';
   /* add hours */
   return time + atof(buf) * 60.;
}

static void bye(void)
{
   exit(0);
}

/* Copy n-1 chars from src[] to dest[].  If src[] is too short, pad dest[]
with leading '0's.  If src[] is too long, not all of src[] will be copied. 
dest[] will always receive n-1 chars followed by null terminator. */
static void cpy0p(char *dest, char *src, int n)
{
   int i, len;

   len = strlen(src);
   for (i = n; --i; )
      *dest++ = (i > len) ? '0' : *src++;
   *dest = '\0';
}

/* Input a double, with default.  Prints name of orbital element
corresponding to elemp, followed by *dp enclosed in brackets.  If no string
is entered, does nothing.  Otherwise, converts entered string to double, adds
value to *dp, and sets "iflag". */

void el_in(int u, struct elemen *elemp, double *dp, int *iflag)
{
   if (u == 2) printf("xndt2o = %11.8lf\n", xndt2o / d1con);

   printf(elemp->elname);

   /* *elemp->conv changes from program units to human units */
   printf(elemp->prnstr, *dp / *elemp->conv);

   if (strlen(s_in("]:", buf))) {
      *dp += atof(buf) * *elemp->conv;
      *iflag = 1;
    }
}

/* Reduces x to range 0 - 2pi */
double fmod2p(double x)
{
   x /= twopi;
   x = (x - dint(x)) * twopi;
   if (x < 0.)
      return x + twopi;
   else
      return x;
}


/* Returns Modified Julian Date (unit = days),
Gregorian calendar.  Any integer is legal for d. */
long int julday(int y, int m, int d)      /* year, month, day */
{
   if (y < 0) {
      printf("Illegal year\n");
      exit(3);
   }

   if (m > 2)   /* after Feb */
      d += leapyr(y);      /* add 1 if leap yr */

   /* Decrement year by 1! */
   y--;

   return (y * 365L + y / 4 - y / 100 + y / 400)
      + calndr[m-1] + d + 1721425 - 2450000;
   /* 1721425 = Julian Date for 1 BC Dec 31 12h, Gregorian calendar */
   /* "Modified Julian Date" by subtracting 2450000 */
}


/* Returns 1 if y is leap year in Gregorian calendar, 0 otherwise. */
static int leapyr(int y)
{
   return((((y % 4) == 0) && ((y % 100) != 0)) || ((y % 400) == 0));
}

/* Prints prompt[] on console, puts typed string in buffer[].  Backspacing
will correct typing mistakes.  Will not allow backspacing past the start of
buffer.  Returns pointer to buffer. */
char *s_in(char *prompt, char *buffer)
{
   printf("%s", prompt);
   if (feof(stdin)) exit(0);
   gets(buffer);
   if (batch) printf("%s\n", buffer);
   return buffer;
}


/* Converts "string" to all upper case, returns "string". */
static char *stoup(char *string)
{
   int c;
   char *cp;

   for (cp = string; (c = *cp) != 0; cp++)
      if (isalpha(c))
         *cp = toupper(c);
   return string;
}


/* Returns Greenwich hour angle of the mean equinox at ep.  If ep is UT1,
returned value is within .001 sec of time compared to the formula in the '85
Astronomical Almanac, for 1990 through 2010.  As a side effect, this function
sets ds50 to days since 1950 Jan 0 0h UT. */
double thetag(double ep)
{
   double ds50;
   double theta;

   ds50 = ep + 16718.0; /* bias 2450000 days (MJD means no 0.5) */
   theta = .27524987 + 1.00273790935 * ds50;   /* revolutions */

   return (theta - dint(theta)) * twopi;
}

/* Tokenizes external char array buf[].  Each ' ' in  buf[] will be replaced
with '\0' & successive members of tokens[] will point to the sub-strings
("tokens") thus created.  End of tokens[] is marked by NULL.  On exit, tokp
points to tokens[0]. */
void tok(void)
{
   char *cptr;
   char c;
   int notok;         /* flag; 1 = not in a token */
   static char *tokens[15];   /* pointers to tokens */

   cptr = buf;
   tokp = tokens;      /* point to 1st element of tokens[] */
   notok = 1;      /* set "not in token" flag true */
   while ((c = *cptr) != '\0') {  /* tokenize the command line */
      if (isspace(c)) {
         notok = 1;
         *cptr = '\0';   /*  replace ' ' with '\0' */
      } else if (notok) {   /* first char of a token */
         *tokp++ = cptr;
         notok = 0;
      }
      cptr++;
   }
   *tokp = NULL;      /* terminate tokens[] */
   tokp = tokens;
}

static void write_elem(char *file)
{

#if 0
/* Satellite desigations */
CLASS char sat_desig[9];
CLASS long sat_ncat;  /* NORAD number */
CLASS int epoch_year;  /* epoch year */

/* satellite's orbital elements. */
CLASS double
   qmo,   /* quasi mean anomaly (omegao + xmo) */
   xmo,   /* mean anomaly */
   xnodeo,   /* right ascension of ascending node */
   omegao,   /* argument of the perigee */
   eo,   /* eccentricity */
   xincl,   /* inclination */
   xno,   /* mean motion, radians/min */
   xndt2o,   /* 1st time derivative of mean motion, or ballistic
      coefficient (depending on ephemeris type) */
   bstar,   /* BSTAR drag term if GP4 theory was used;
      otherwise, radiation pressure coefficient */
   epoch_day;   /* epoch day of elements */
#endif

/* USA 101 */
/* 1 23030U 94 17  A 94221.97749609 0.00000500  00000-0  38850-4 0    04 */
/* 2 23030 105.0307 124.1470 0007385 138.4686 221.5313 15.01273715    03 */

   int year;
   double xno2;
   char xno2_string[11];
   char bstar_string[13];
   char bstar_fract[7];
   char bstar_exp[3];
   char eo_string[11];

/*
   strncpy(sat_desig, &input_line[9], 8);
*/

   xno2 = xndt2o / d1con;
   sprintf(xno2_string, "%10.8lf", fabs(xno2));
   if (xno2_string[0] == '0') xno2_string[0] = ' ';
   if (xno2 < 0.0) xno2_string[0] = '-';

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

/* printf(" epoch_year %d epoch_day %11.8f\n", epoch_year, epoch_day); */

   year = epoch_year - 1900;
   if (year >= 100) year -= 100;
   sprintf(line1, "1 %05ldU %8s %02d%012.8lf %10s  00000-0 %6s%2s 0    00",
           sat_ncat, sat_desig, year, epoch_day, xno2_string,
            bstar_fract, bstar_exp);
   line1[68] = cksum(line1) + '0';

   sprintf(eo_string, "%10.7lf", eo);
   xnodeo = fmod2p(xnodeo);
   omegao = fmod2p(omegao);
   xmo = fmod2p(qmo - omegao);
   sprintf(line2, "2 %05ld %8.4lf %8.4lf %7s %8.4lf %8.4lf %11.8lf    00",
           sat_ncat, xincl * ra2de, xnodeo * ra2de, &eo_string[3],
           omegao * ra2de, xmo * ra2de, xno / twopi * 1440.0);
   line2[68] = cksum(line2) + '0';

   printf("\n%s\n", sat_name);
   printf("%s\n", line1);
   printf("%s\n", line2);

   if(strlen(file) > 0)
   {
      FILE *fp;
      if(out) fp = fopen(file, "w");
      else
      {
         fp = fopen(file, "a");
         fprintf(fp, "\n");
      }
      fprintf(fp, "%s\n", sat_name);
      fprintf(fp, "%s\n", line1);
      fprintf(fp, "%s\n", line2);
      fclose(fp);
   }
}

