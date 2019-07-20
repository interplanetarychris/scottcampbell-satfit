/* getel.c
by Paul S. Hirose, 1990 Sep 27
Handles orbital element files and reads the elements.
*/

#include "elcor.h"   /* global header */

/*######################## LOCAL FUNCTIONS ########################*/

/* convert string to angle in radians */
static double angle(char *string);
/* checksums a line of data */
int cksum(char *line);
/* converts NORAD Epoch string to epoch_year and epoch_day */
static void epjd(char *buf);
/* replace spaces in buffer with nulls */
static void nullsp(char *buf);
/* returns value of scientific notation #'s */
static double sci(char *string);
/* skip leading nulls */
static char *sknull(char *cptr);

/*############################## CODE ##############################*/

/* Returns the radian value of the angle (degrees) in "string".
Leading nulls on "string" will be skipped, so it had better contain
at least one digit. */
static double angle(char *string)
{
   return atof(sknull(string)) * de2ra;
}


/* Does a checksum modulo 10 on the given line.
   Digits = their value, '-' = 1, all other chars = 0.
   Returns 0 if ok. */
int cksum(char *line)
{
   int tot = 0;   /* checksum accumulator */
   int count = 68;      /* length of data line */

   /* accumulate checksum from all but the last char in the line,
   which is the desired checksum */
   do {
      int c;

      c = *line++;
      if (isdigit(c))
         tot += c - '0';
      else if (c == '-')
         tot++;
      /* all other chars = 0 */
   } while (--count);

   /* Convert accumulated # to mod 10, subtract desired sum */
   return((tot % 10) - (*line - '0'));
}

/* Given pointer to the orbital elements Epoch field, sets epoch_year
and epoch_day in global area.  Format of Epoch field:  YYDDD.dddddddd
*/
static void epjd(char *buf)
{
   epoch_day = atof(sknull(buf + 2));

   buf[2] = '\0';         /* day 100s digit = null */
   epoch_year = atoi(sknull(buf));
   epoch_year += (epoch_year >= 57) ? 1900 : 2000;
}

/* Sets the global orbital element variables from file fp.  Initializes
precession.  Prints warning if checksum error in
file, but will still attempt to use elements. */
int getel(FILE *fp)
{
   int len;
#if 0
   double ao, a1, delo, del1, temp;
#endif
   int check;

   fgets(buf, 80, fp);   /* read the name */
   buf[80] = '\0';
   strcpy(sat_name, buf);
   len = strlen(sat_name);
   if (sat_name[len-1] == '\n')
      sat_name[len-1] = '\0';

   /* Get first data line from file */
   if (fgets(line1, 80, fp) == NULL)
      return -1;   /* incorrect format */
   if (line1[0] != '1')
      return -1;   /* incorrect format */
   check = cksum(line1);   /* Line 1 checksum */

   memcpy(sat_desig, &line1[9], 8);
   sat_desig[8] = '\0';

   nullsp(line1);

   sat_ncat = atol(sknull(&line1[2]));

   if (line1[35])      /* field isn't blank */
      /* multiply xndt20 by 2pi / (1440^2) */
      xndt2o = atof(sknull(line1 + 33)) * d1con;
   else
      xndt2o = 0.;
   bstar = sci(line1 + 53);
   epjd(line1 + 18);


   /* get second data line */
   if (fgets(line2, 80, fp) == NULL)
      return -1;   /* incorrect format */
   if (line2[0] != '2')
      return -1;   /* incorrect format */
   check += cksum(line2);   /* add Line 2 checksum */

   nullsp(line2);      /* turn spaces to nulls */

   xmo = angle(line2 + 43);
   xnodeo = angle(line2 + 17);
   omegao = angle(line2 + 34);
   qmo = omegao + xmo;      /* quasi mean anomaly */
   xincl = angle(line2 + 8);
   line2[25] = '.';      /* add decimal point to eccentricity */
   eo = atof(line2 + 25);

   /* Make sure mean motion is null-terminated, since rev. no.
   may immediately follow. */
   buf[0] = line2[63];
   line2[63] = '\0';
   /* multiply xno by 2pi / 1440 */
   xno = atof(sknull(line2 + 52)) * nocon;
   line2[63] = buf[0];

   if (check)
      printf("CAUTION:  FAILED CHECKSUM\n");

#if 0   /* leave this stuff out for now */

   /* Compute period.  Object is deep space if <= 6.4 revs/day */
   a1 = POW(xke / xno, tothrd);
   temp = cos(xincl);
   temp = 1.5 * ck2 * (3. * temp * temp - 1.) /
     POW(1. - eo * eo, 1.5);
   del1 = temp / (a1 * a1);
   ao = a1 * (1. - del1 * (.5 * tothrd + del1 *
     (1. + 134. / 81. * del1)));
   delo = temp / (ao * ao);

   if (xno / (1. + delo) <= .0279253)
      return 1;   /* deep space */
   else
#endif
      return 0;   /* near earth */
}


/* replaces each occurrence of ' ' in string "buf" with a null */
static void nullsp(char *buf)
{
   char c;

   while ((c = *buf) != '\0') {
      if (c == ' ')
         *buf = '\0';
      buf++;
    }
}


/* Converts quasi scientific notation to double.  Format:  SnSn where S is
either '-' or '\0' and n represents a sequence of 1 or more digits.  An
implied decimal point exists to the left of the left n.  Total length 27
chars max. */
static double sci(char *string)
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


/* Skip nulls.  Returns pointer to first non-null char in "cptr". */
static char *sknull(char *cptr)
{
   while (!*cptr)
      ++cptr;

   return cptr;
}
